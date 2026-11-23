// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of one output grouped by input record, and the input lines that tell two
//! records under one identifier apart.

use super::entry::{replaced_reference, uploaded_alleles, Entry, Fields, Row};
use super::shared::dash;
use super::Refusal;
use crate::pytext;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

/// The largest coordinate vep-rs holds; a POS past it fails its line in every input format.
pub const VEP_RS_COORDINATE_CEILING: i128 = (1 << 62) - 1;

/// The input record an entry belongs to: the Location with the `Uploaded_variation`, which
/// tells two records at one Location apart (two records sharing both are read as one), and
/// the first segment of the key alone for an entry without a Location.
pub fn record_label(e: &Entry) -> String {
    let label = e.key.split('|').next().unwrap_or("");
    if e.location.is_none() {
        return label.to_string();
    }
    match dash(e.fields.get_or("Uploaded_variation", "-")) {
        "-" => label.to_string(),
        uploaded => format!("{label} {uploaded}"),
    }
}

/// The fields Perl's list-context `split` keeps: the trailing empty fields dropped and the
/// middle ones kept, so `G,,T,,` is `G`, the empty field, `T`.
pub fn perl_fields(mut fields: Vec<&str>) -> Vec<&str> {
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

/// Perl's list-context `split` of `text` on the separator `sep`, as Ensembl VEP splits the
/// ALT column (BaseVCF4.pm:489): `<DEL>,` is the one field `<DEL>` and the empty string has
/// no field.
pub fn perl_split<'a>(text: &'a str, sep: &str) -> Vec<&'a str> {
    perl_fields(text.split(sep).collect())
}

/// Whether `text` is ASCII digits with an optional fraction and exponent and an optional
/// `+` before them, the spellings of a number that Rust's `f64` parser accepts.
fn pos_number(text: &str) -> bool {
    let unsigned = text.strip_prefix('+').unwrap_or(text);
    let (mantissa, exponent) = match unsigned.find(['e', 'E']) {
        Some(at) => (&unsigned[..at], Some(&unsigned[at + 1..])),
        None => (unsigned, None),
    };
    let digits = |part: &str| part.bytes().all(|b| b.is_ascii_digit());
    let mantissa_read = match mantissa.split_once('.') {
        Some((whole, fraction)) => {
            digits(whole) && digits(fraction) && whole.len() + fraction.len() > 0
        }
        None => pytext::is_ascii_digits(mantissa),
    };
    mantissa_read
        && exponent.is_none_or(|exponent| {
            pytext::is_ascii_digits(exponent.strip_prefix(['+', '-']).unwrap_or(exponent))
        })
}

/// The POS column as vep-rs's parsers read it: an unpadded decimal integer up to
/// `VEP_RS_COORDINATE_CEILING` exactly, or another spelling of an exact non-negative integer
/// below the ceiling that Rust's `f64` parser accepts, blanks around it ignored (`3e7`,
/// `30425650.0`, `30425650.`, `+30425650`), read as the double it names; `None` where vep-rs
/// fails the line (`.`, `abc`, `-5`, `0x1D04A1A`, `30_425_650`, digits of another script,
/// `inf`, `1e400`, `30425650.5`).
pub fn vep_rs_pos(text: &str) -> Option<i128> {
    let stripped = text.trim_matches([' ', '\t', '\n', '\r', '\x0c']);
    if !pos_number(stripped) {
        return None;
    }
    if pytext::isdigit(text.strip_prefix('+').unwrap_or(text)) {
        let value = pytext::int(text)?;
        return (value <= VEP_RS_COORDINATE_CEILING).then_some(value);
    }
    // The ceiling rounds up to 2^62 as a double, so the bound is strict.
    let number: f64 = stripped.parse().ok()?;
    (number.is_finite() && number.fract() == 0.0 && number < 2f64.powi(62))
        .then_some(number as i128)
}

/// Whether the input line `columns` (`CHROM`, `POS`, `ID`, `REF`, `ALT`) writes the default
/// or tab row `e`: it stands on the row's chromosome, a `chr` prefix set aside, and one of
/// its alternates is its reference allele with the row's allele in place of the Location's
/// bases, or between the Location's two positions (`replaced_reference`). Each alternate is
/// asked by itself, so a line two of whose alternates give the row in two ways writes it too.
pub fn line_writes_row(columns: &[String], e: &Entry) -> Result<bool, Refusal> {
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    let Some(position) = vep_rs_pos(&columns[1]) else {
        return Ok(false);
    };
    let without_chr = |name: &str| name.strip_prefix("chr").unwrap_or(name).to_string();
    if without_chr(&columns[0]) != without_chr(pytext::partition(location, ":").0) {
        return Ok(false);
    }
    for alternate in perl_split(&columns[4], ",") {
        let read = replaced_reference(
            position,
            &columns[3],
            &[alternate],
            Some(location),
            &e.allele,
        )?;
        if read.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The rows of one label under the number of the input line that writes each.
type RowsByLine = HashMap<usize, Vec<Row>>;

/// The entries of one output grouped by input record (`record_label`): the row group a
/// class reads whose shape spans a record's rows. The input lines are the input VCF's data
/// lines, read for the rows that name their record by an identifier. A line carries the
/// first `;`-joined element of its ID column, the name both engines print for its rows, and
/// no name where that element is a dot, empty or `0`.
pub struct RecordRows {
    by_record: HashMap<String, Vec<Row>>,
    /// The named input lines, each with the identifier as its third column.
    lines: Vec<Vec<String>>,
    lines_by_name: HashMap<String, Vec<usize>>,
    /// The entries the rows were built from, by address.
    members: HashSet<*const Entry>,
    /// The line of each member that several lines were asked about.
    line_by_entry: RefCell<HashMap<*const Entry, Option<usize>>>,
    /// Per label, its entries under the line that writes each, in the label's order.
    entries_by_line: RefCell<HashMap<String, Rc<RowsByLine>>>,
}

impl RecordRows {
    pub fn new<S: AsRef<str>>(entries: &[Row], input_lines: &[S]) -> RecordRows {
        let mut by_record: HashMap<String, Vec<Row>> = HashMap::new();
        for e in entries {
            by_record
                .entry(record_label(e))
                .or_default()
                .push(Rc::clone(e));
        }
        let mut lines = Vec::new();
        let mut lines_by_name: HashMap<String, Vec<usize>> = HashMap::new();
        for line in input_lines {
            let mut columns: Vec<String> = line.as_ref().split('\t').map(str::to_string).collect();
            if columns.len() < 5 {
                continue;
            }
            let name = columns[2].split(';').next().unwrap_or("").to_string();
            if matches!(name.as_str(), "." | "" | "0") {
                continue;
            }
            columns[2] = name.clone();
            lines_by_name.entry(name).or_default().push(lines.len());
            lines.push(columns);
        }
        RecordRows {
            by_record,
            lines,
            lines_by_name,
            members: entries.iter().map(Rc::as_ptr).collect(),
            line_by_entry: RefCell::new(HashMap::new()),
            entries_by_line: RefCell::new(HashMap::new()),
        }
    }

    /// The entries under `label`, in the order they were given.
    pub fn labelled(&self, label: &str) -> &[Row] {
        self.by_record.get(label).map_or(&[], Vec::as_slice)
    }

    /// The number of the input line of the record a row under an identifier belongs to: the
    /// one line that carries the identifier, or of several the one that writes the row at
    /// its Location with its allele (`line_writes_row`). `None` for an identifier no line
    /// carries, and where none or more than one of several lines writes the row.
    fn line_number(&self, e: &Row) -> Result<Option<usize>, Refusal> {
        let name = e.fields.get_or("Uploaded_variation", "");
        let lines = self.lines_by_name.get(name).map_or(&[][..], Vec::as_slice);
        if lines.len() < 2 {
            return Ok(lines.first().copied());
        }
        let address = Rc::as_ptr(e);
        if let Some(asked) = self.line_by_entry.borrow().get(&address) {
            return Ok(*asked);
        }
        let mut writing = Vec::new();
        for &number in lines {
            if line_writes_row(&self.lines[number], e)? {
                writing.push(number);
            }
        }
        let line = if writing.len() == 1 {
            Some(writing[0])
        } else {
            None
        };
        if self.members.contains(&address) {
            self.line_by_entry.borrow_mut().insert(address, line);
        }
        Ok(line)
    }

    /// The input line of the record a row under an identifier belongs to (`CHROM`, `POS`,
    /// the identifier, `REF`, `ALT` and the columns after them), as `line_number` finds it.
    pub fn line_of(&self, e: &Row) -> Result<Option<&[String]>, Refusal> {
        Ok(self
            .line_number(e)?
            .map(|number| self.lines[number].as_slice()))
    }

    /// The entries of `label` whose line is the line `number`, in the label's order.
    fn written_by(&self, label: &str, number: usize) -> Result<Vec<Row>, Refusal> {
        let held = self.entries_by_line.borrow().get(label).cloned();
        let by_line = match held {
            Some(by_line) => by_line,
            None => {
                let mut by_line = RowsByLine::new();
                for row in self.labelled(label) {
                    if let Some(written) = self.line_number(row)? {
                        by_line.entry(written).or_default().push(Rc::clone(row));
                    }
                }
                let by_line = Rc::new(by_line);
                self.entries_by_line
                    .borrow_mut()
                    .insert(label.to_string(), Rc::clone(&by_line));
                by_line
            }
        };
        Ok(by_line.get(&number).cloned().unwrap_or_default())
    }

    /// Whether `e`'s label names one input record: an entry without a Location, a row whose
    /// identifier-less `Uploaded_variation` spells the record's alleles, or one whose
    /// identifier names one line of the input (`line_of`). Under an identifier two records
    /// at one Location share a label, so without the input their rows do not tell a
    /// multi-allelic record from two one-allele records.
    pub fn one_record(&self, e: &Row) -> Result<bool, Refusal> {
        if e.location.is_none() {
            return Ok(true);
        }
        if uploaded_alleles(e.fields.get_or("Uploaded_variation", "")).is_some() {
            return Ok(true);
        }
        Ok(self.line_number(e)?.is_some())
    }

    /// The entries of `e`'s input record, or `None` where the label may name two records
    /// (`one_record`). For a row under an identifier an entry built from its input line
    /// (`CHROM`, `POS`, `REF`, `ALT`) follows the rows, a new entry at every call. Where
    /// several lines carry the identifier, the rows are those of the label that the same
    /// line writes.
    pub fn of(&self, e: &Row) -> Result<Option<Vec<Row>>, Refusal> {
        if !self.one_record(e)? {
            return Ok(None);
        }
        let label = record_label(e);
        let line = match e.location {
            Some(_) => self.line_number(e)?,
            None => None,
        };
        let Some(number) = line else {
            return Ok(Some(self.labelled(&label).to_vec()));
        };
        let line = &self.lines[number];
        let mut rows = if self.lines_by_name[&line[2]].len() > 1 {
            self.written_by(&label, number)?
        } else {
            self.labelled(&label).to_vec()
        };
        let fields: Fields = [
            ("Uploaded_variation", &line[2]),
            ("CHROM", &line[0]),
            ("POS", &line[1]),
            ("REF", &line[3]),
            ("ALT", &line[4]),
        ]
        .into_iter()
        .collect();
        rows.push(Rc::new(Entry {
            key: format!("{}|record", e.location.as_deref().unwrap_or("")),
            location: e.location.clone(),
            fields,
            ..Entry::default()
        }));
        Ok(Some(rows))
    }
}
