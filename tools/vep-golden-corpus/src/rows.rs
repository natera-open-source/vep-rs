// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of a default-format output and the lines of an input VCF as the subcommands
//! read them, and the consequence set a row is keyed by.

pub mod blocks;
pub mod text;

pub use text::{not_readable_text, ByteSource, ReadError, TextLines};

use crate::classes::{
    record_class_function, row_family_classes, variant_class, FamilySpec, RECORD_CLASS_RESERVOIR,
};
use crate::error::{failure, Failure};
use crate::names::{length, names_of_record, record_pos};
use crate::pyfmt::py_strip;
use crate::pytext::{path_str, py_int, py_repr_str};
use memchr::memmem::Finder;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

/// `opener(path)(path, "rt", encoding="utf-8", errors="replace")`: the lines of a text
/// file, read through `gzip.open` for a name ending in `.gz` or `.bgz` (BGZF is gzip) and
/// through `open` otherwise.
pub fn opener(path: &str) -> Result<TextLines, ReadError> {
    TextLines::open(path)
}

/// A Consequence column as a key holds it: its terms without repeats, sorted by code point
/// and joined by `,`; blanks around the column and empty terms are dropped.
pub fn normalize_consequence_set(raw: &str) -> String {
    let terms: BTreeSet<&str> = py_strip(raw)
        .split(',')
        .filter(|term| !term.is_empty())
        .collect();
    terms.into_iter().collect::<Vec<&str>>().join(",")
}

/// `normalize_consequence_set(raw)`, the column itself where it is written that way.
fn consequence_set(raw: &str) -> Cow<'_, str> {
    let mut previous = "";
    let rising = raw.split(',').all(|term| {
        let after = previous < term;
        previous = term;
        after
    });
    if raw.is_empty() || (rising && py_strip(raw).len() == raw.len()) {
        Cow::Borrowed(raw)
    } else {
        Cow::Owned(normalize_consequence_set(raw))
    }
}

/// The fields of an Extra column: a key holds the value it was given last and stands where
/// it was given first, as in a `dict`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Extra<'a> {
    fields: Vec<(&'a str, &'a str)>,
}

impl<'a> Extra<'a> {
    /// The value under `key`, `None` for a key the column lacks: `extra.get(key)`.
    pub fn get(&self, key: &str) -> Option<&'a str> {
        self.fields
            .iter()
            .find(|(name, _)| *name == key)
            .map(|(_, value)| *value)
    }

    /// The keys and values in the order of the keys' first places: `extra.items()`.
    pub fn iter(&self) -> impl Iterator<Item = (&'a str, &'a str)> + '_ {
        self.fields.iter().copied()
    }

    pub fn len(&self) -> usize {
        self.fields.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }
}

/// The fields of an Extra column, `KEY=VALUE` items joined by `;`: the value of an item is
/// what follows its first `=`, an item without `=` is a key with an empty value, an empty
/// item is none, and an empty column or `-` has no field.
pub fn parse_extra(raw: &str) -> Extra<'_> {
    let mut fields: Vec<(&str, &str)> = Vec::new();
    if raw.is_empty() || raw == "-" {
        return Extra { fields };
    }
    for part in raw.split(';') {
        let (key, value) = match part.split_once('=') {
            Some(pair) => pair,
            None if part.is_empty() => continue,
            None => (part, ""),
        };
        match fields.iter_mut().find(|(name, _)| *name == key) {
            Some(field) => field.1 = value,
            None => fields.push((key, value)),
        }
    }
    Extra { fields }
}

/// `parse_extra(raw).get(key)` read without the other fields: the value of the last item
/// whose key it is. The key is one whose end is not also its start, so that no two places
/// it could stand at overlap.
fn extra_value<'a>(raw: &'a str, key: &Finder<'_>) -> Option<&'a str> {
    let bytes = raw.as_bytes();
    let mut found = None;
    for at in key.find_iter(bytes) {
        if at > 0 && bytes[at - 1] != b';' {
            continue;
        }
        let after = at + key.needle().len();
        match bytes.get(after) {
            None | Some(b';') => found = Some(""),
            Some(b'=') => found = raw[after + 1..].split(';').next(),
            Some(_) => {}
        }
    }
    found
}

/// A dataset of a command line as `--dataset` names it: its id, its reference output and its
/// input VCF, the two paths as `str(pathlib.Path(...))` writes them.
pub type Dataset = (String, String, String);

/// `sid=reference.txt=input.vcf[.gz]` as a dataset; an argument of another count of `=`
/// is refused.
pub fn parse_dataset_arg(arg: &str) -> Result<Dataset, Failure> {
    let parts: Vec<&str> = arg.split('=').collect();
    match parts[..] {
        [sid, reference, input] => Ok((sid.to_string(), path_str(reference), path_str(input))),
        _ => Err(failure(format!(
            "--dataset wants sid=reference=input, got {}",
            py_repr_str(arg)
        ))),
    }
}

/// The header lines and the records of a VCF (`vcf_lines`).
pub struct VcfLines {
    path: String,
    columns: Option<usize>,
    lines: Option<TextLines>,
    failure: Option<Failure>,
}

/// Every header line (`true`, every column) and every record of at least eight fields
/// (`false`) of a VCF, split on tabs. A header line is the `#CHROM` line, or any other line
/// that starts with one `#`; `##` lines are skipped. With `columns` a record is split into
/// that many fields and one more holding the rest of the line, so a source of thousands of
/// samples is not split beyond the columns read. A file that cannot be opened or read (a
/// `.gz` that is not gzip among them) is refused with its path, in place of the line the
/// reading stopped at.
pub fn vcf_lines(path: &str, columns: Option<usize>) -> VcfLines {
    let (lines, failure) = match opener(path) {
        Ok(lines) => (Some(lines), None),
        Err(error) => (None, Some(not_readable_text(path, &error))),
    };
    VcfLines {
        path: path.to_string(),
        columns,
        lines,
        failure,
    }
}

impl VcfLines {
    /// The next header line or record as the reader's line holds it: whether it is a
    /// header, and its columns. The iterator gives the same lines with columns of their own.
    pub fn next_columns(&mut self) -> Option<Result<(bool, Vec<&str>), Failure>> {
        if let Some(failure) = self.failure.take() {
            return Some(Err(failure));
        }
        loop {
            let lines = self.lines.as_mut()?;
            match lines.advance() {
                Ok(true) => {}
                Ok(false) => {
                    self.lines = None;
                    return None;
                }
                Err(error) => {
                    self.lines = None;
                    return Some(Err(not_readable_text(&self.path, &error)));
                }
            }
            let line = lines.line();
            if line.starts_with("##") {
                continue;
            }
            let tabs = memchr::memchr_iter(b'\t', line.as_bytes()).take(7).count();
            if line.starts_with('#') || tabs.min(self.columns.unwrap_or(7)) == 7 {
                break;
            }
        }
        let line = self.lines.as_ref()?.line();
        let columns = match self.columns {
            Some(columns) if !line.starts_with('#') => {
                line.splitn(columns.saturating_add(1), '\t').collect()
            }
            _ => line.split('\t').collect(),
        };
        Some(Ok((line.starts_with('#'), columns)))
    }
}

impl Iterator for VcfLines {
    type Item = Result<(bool, Vec<String>), Failure>;

    fn next(&mut self) -> Option<Self::Item> {
        let (is_header, columns) = match self.next_columns()? {
            Ok(line) => line,
            Err(failure) => return Some(Err(failure)),
        };
        let columns = columns.into_iter().map(str::to_string).collect();
        Some(Ok((is_header, columns)))
    }
}

/// The source `kept_samples` names when the samples to agree with are an earlier input's.
pub const EARLIER_INPUT: &str = "an earlier input";

/// The first `keep` sample names of a `#CHROM` header line, refused when the line lacks
/// `FORMAT` or that many samples, or when `samples` (the names `source` keeps, when it
/// keeps any) names others.
pub fn kept_samples<S: AsRef<str>>(
    header: &[S],
    keep: usize,
    path: &str,
    samples: &[String],
    source: &str,
) -> Result<Vec<String>, Failure> {
    if header.len() < 9 || header[8].as_ref() != "FORMAT" {
        return Err(failure(format!(
            "{path} has no FORMAT column, so no sample columns"
        )));
    }
    if header.len() < 9 + keep {
        return Err(failure(format!(
            "{path} has {} sample columns, {keep} are kept",
            header.len() - 9
        )));
    }
    let kept: Vec<String> = header[9..9 + keep]
        .iter()
        .map(|name| name.as_ref().to_string())
        .collect();
    if !samples.is_empty() && kept.as_slice() != samples {
        return Err(failure(format!(
            "{path} keeps samples {} where {source} keeps {}: every input must keep the same samples",
            kept.join(","),
            samples.join(",")
        )));
    }
    Ok(kept)
}

/// One row of a default-format output.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Row {
    pub dataset: String,
    pub uploaded_variation: String,
    pub location: String,
    pub allele: String,
    pub feature: String,
    pub feature_type: String,
    /// The Consequence column as `normalize_consequence_set` writes it.
    pub consequence_set: String,
    /// The Extra field `STRAND`, empty when the row has none.
    pub strand: String,
    /// The Extra field `FLAGS`, empty when the row has none.
    pub flags: String,
    /// The columns cDNA_position, CDS_position and Protein_position.
    pub positions: [String; 3],
    /// The row's class under each spec of `--family-columns`.
    pub family: Vec<String>,
    /// Whether the input record holds the row's allele as a substitution: what tells a
    /// two-base substitution from a two-base insertion.
    pub substitution: bool,
}

/// One row of a default-format output as the reader's line holds it: the members of `Row`
/// and, on a line that has all twelve columns up to Codons, the columns Protein_position,
/// Amino_acids and Codons, borrowed until the reader moves on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RowView<'a> {
    pub dataset: &'a str,
    pub uploaded_variation: &'a str,
    pub location: &'a str,
    pub allele: &'a str,
    pub feature: &'a str,
    pub feature_type: &'a str,
    pub consequence_set: Cow<'a, str>,
    pub strand: &'a str,
    pub flags: &'a str,
    pub positions: [&'a str; 3],
    /// The Protein_position, Amino_acids and Codons columns of a line that has all twelve
    /// columns up to Codons.
    pub cells: Option<[&'a str; 3]>,
    pub family: Vec<String>,
    pub substitution: bool,
}

/// `;`-joined, the family classes and the variant class.
fn family_stratum(family: &[String], klass: &str) -> String {
    let mut parts: Vec<&str> = family.iter().map(String::as_str).collect();
    parts.push(klass);
    parts.join(";")
}

/// How much a row is preferred as an exemplar: 4 for a row with FLAGS, 3 for a position
/// that holds `?`, 2 for a position that is a range, summed.
fn preference(flags: &str, positions: [&str; 3]) -> i64 {
    let mut score = 0;
    if !flags.is_empty() {
        score += 4;
    }
    if positions.iter().any(|position| position.contains('?')) {
        score += 3;
    }
    if positions
        .iter()
        .any(|position| position.contains('-') && *position != "-")
    {
        score += 2;
    }
    score
}

/// The genomic width of a Location in bases: 0 for a single position and for coordinates
/// `int()` does not read.
fn span(location: &str) -> i64 {
    let coordinates = location.split_once(':').map_or(location, |(_, tail)| tail);
    let Some((start, end)) = coordinates.split_once('-') else {
        return 0;
    };
    match (py_int(start), py_int(end)) {
        (Some(start), Some(end)) => {
            let width = (i128::from(end) - i128::from(start)).abs();
            i64::try_from(width).unwrap_or(i64::MAX)
        }
        _ => 0,
    }
}

impl Row {
    /// What the rows of one input record share: Uploaded_variation and Location.
    pub fn record_key(&self) -> (&str, &str) {
        (&self.uploaded_variation, &self.location)
    }

    /// The row's variant class (`variant_class`).
    pub fn klass(&self) -> String {
        variant_class(&self.location, &self.allele, self.substitution)
    }

    /// The family classes and the variant class, `;`-joined: the key `--family-columns`
    /// covers.
    pub fn family_stratum(&self) -> String {
        family_stratum(&self.family, &self.klass())
    }

    /// How much the row is preferred as an exemplar: 4 for a row with FLAGS, 3 for a
    /// position that holds `?`, 2 for a position that is a range, summed.
    pub fn preference(&self) -> i64 {
        let [cdna, cds, protein] = &self.positions;
        preference(&self.flags, [cdna, cds, protein])
    }

    /// The genomic width of the Location in bases, 0 for a single position.
    pub fn span(&self) -> i64 {
        span(&self.location)
    }
}

impl RowView<'_> {
    /// The row with members of its own.
    pub fn to_row(&self) -> Row {
        Row {
            dataset: self.dataset.to_string(),
            uploaded_variation: self.uploaded_variation.to_string(),
            location: self.location.to_string(),
            allele: self.allele.to_string(),
            feature: self.feature.to_string(),
            feature_type: self.feature_type.to_string(),
            consequence_set: self.consequence_set.to_string(),
            strand: self.strand.to_string(),
            flags: self.flags.to_string(),
            positions: self.positions.map(str::to_string),
            family: self.family.clone(),
            substitution: self.substitution,
        }
    }

    /// `Row::record_key`.
    pub fn record_key(&self) -> (&str, &str) {
        (self.uploaded_variation, self.location)
    }

    /// `Row::klass`.
    pub fn klass(&self) -> String {
        variant_class(self.location, self.allele, self.substitution)
    }

    /// `Row::family_stratum`.
    pub fn family_stratum(&self) -> String {
        family_stratum(&self.family, &self.klass())
    }

    /// `Row::preference`.
    pub fn preference(&self) -> i64 {
        preference(self.flags, self.positions)
    }

    /// `Row::span`.
    pub fn span(&self) -> i64 {
        span(self.location)
    }
}

/// Values by chromosome and Uploaded_variation name: what a `dict` keyed by
/// `(chromosome, name)` holds, read without building a key.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ByRecordName<T> {
    by_chromosome: HashMap<String, HashMap<String, Vec<T>>>,
}

impl<T> Default for ByRecordName<T> {
    fn default() -> Self {
        ByRecordName {
            by_chromosome: HashMap::new(),
        }
    }
}

impl<T> ByRecordName<T> {
    /// The values under `(chromosome, name)`, `None` for a key that holds none:
    /// `d.get((chromosome, name))`.
    pub fn get(&self, chromosome: &str, name: &str) -> Option<&[T]> {
        let names = self.by_chromosome.get(chromosome)?;
        names.get(name).map(Vec::as_slice)
    }

    /// Adds values under `(chromosome, name)`, after the ones it holds:
    /// `d.setdefault((chromosome, name), []).extend(values)`.
    pub fn extend(&mut self, chromosome: &str, name: &str, values: impl IntoIterator<Item = T>) {
        let names = match self.by_chromosome.get_mut(chromosome) {
            Some(names) => names,
            None => self
                .by_chromosome
                .entry(chromosome.to_string())
                .or_default(),
        };
        match names.get_mut(name) {
            Some(held) => held.extend(values),
            None => {
                names.insert(name.to_string(), values.into_iter().collect());
            }
        }
    }

    /// The count of keys.
    pub fn len(&self) -> usize {
        self.by_chromosome.values().map(HashMap::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// The (REF, ALT) pairs of the records of an input VCF (`input_alleles`).
pub type InputAlleles = ByRecordName<(String, String)>;

/// The POS and bases of the substitution ALTs of an input VCF (`substitution_alts`).
pub type Substitutions = ByRecordName<(i64, String)>;

/// Every (REF, ALT) pair of an input VCF's records by (chromosome, name), under every
/// Uploaded_variation name the record can carry (`record_names`).
pub fn input_alleles(path: &str) -> Result<InputAlleles, Failure> {
    let mut out = InputAlleles::default();
    let mut lines = vcf_lines(path, None);
    while let Some(line) = lines.next_columns() {
        let (is_header, f) = line?;
        if is_header {
            continue;
        }
        let pairs = || {
            f[4].split(',')
                .map(|alt| (f[3].to_string(), alt.to_string()))
        };
        for name in names_of_record(f[0], f[1], f[2], f[3], f[4])? {
            out.extend(f[0], &name, pairs());
        }
    }
    Ok(out)
}

/// The POS and bases of every ALT of an input VCF that is as long as its record's REF of
/// two bases or more, a substitution of that many bases, by (chromosome, name), under every
/// Uploaded_variation name the record can carry (`record_names`). The ALT is ASCII letters.
pub fn substitution_alts(path: &str) -> Result<Substitutions, Failure> {
    let mut out = Substitutions::default();
    let mut lines = vcf_lines(path, Some(8));
    while let Some(line) = lines.next_columns() {
        let (is_header, f) = line?;
        if is_header || length(f[3]) < 2 {
            continue;
        }
        let width = length(f[3]);
        let letters = |alt: &&str| !alt.is_empty() && alt.bytes().all(|b| b.is_ascii_alphabetic());
        let alts: Vec<&str> = f[4]
            .split(',')
            .filter(|alt| length(alt) == width)
            .filter(letters)
            .collect();
        if alts.is_empty() {
            continue;
        }
        let pos = record_pos(f[0], f[1], f[2])?;
        for name in names_of_record(f[0], f[1], f[2], f[3], f[4])? {
            out.extend(f[0], &name, alts.iter().map(|alt| (pos, alt.to_string())));
        }
    }
    Ok(out)
}

/// Whether a substitution of `alts` (`substitution_alts`) holds `allele` at the positions
/// of `span`, a row's Location after its chromosome: the row is then that substitution's,
/// whole or as Ensembl VEP prints it once the bases its alleles share are dropped (the
/// anchor base of a multi-allelic record with an indel, either end under `--minimal`). The
/// record's position is part of the test, so a row never reads a record that only shares
/// its name. A span `int()` does not read holds no allele.
pub fn holds_allele(span: &str, allele: &str, alts: &[(i64, String)]) -> bool {
    let (start, end) = span.split_once('-').unwrap_or((span, ""));
    let end = if end.is_empty() { start } else { end };
    let (Some(a), Some(b)) = (py_int(start), py_int(end)) else {
        return false;
    };
    let (a, b) = (i128::from(a), i128::from(b));
    let Ok(width) = usize::try_from(b - a + 1) else {
        return false;
    };
    if length(allele) != width {
        return false;
    }
    let wanted = allele.to_uppercase();
    alts.iter().any(|(pos, alt)| {
        let Ok(skip) = usize::try_from(a - i128::from(*pos)) else {
            return false;
        };
        let held: String = alt.chars().skip(skip).take(width).collect();
        held.to_uppercase() == wanted
    })
}

/// The Extra keys a row's members are read from.
struct RowKeys {
    strand: Finder<'static>,
    flags: Finder<'static>,
}

/// The row a line of at least seven columns holds; a column the line lacks is `-`.
fn row_view<'a>(
    line: &'a str,
    dataset: &'a str,
    keys: &RowKeys,
    specs: &[FamilySpec],
    substitutions: Option<&Substitutions>,
) -> RowView<'a> {
    let mut f = ["-"; 14];
    let mut held = 0;
    for (column, text) in f.iter_mut().zip(line.split('\t')) {
        *column = text;
        held += 1;
    }
    let (chromosome, positions) = f[1].split_once(':').unwrap_or((f[1], ""));
    let alts = substitutions.and_then(|held| held.get(chromosome, f[0]));
    RowView {
        dataset,
        uploaded_variation: f[0],
        location: f[1],
        allele: f[2],
        feature: f[4],
        feature_type: f[5],
        consequence_set: consequence_set(f[6]),
        strand: extra_value(f[13], &keys.strand).unwrap_or(""),
        flags: extra_value(f[13], &keys.flags).unwrap_or(""),
        positions: [f[7], f[8], f[9]],
        cells: (held >= 12).then_some([f[9], f[10], f[11]]),
        family: if specs.is_empty() {
            Vec::new()
        } else {
            row_family_classes(f[13], specs, f[1], f[2])
        },
        substitution: alts.is_some_and(|alts| holds_allele(positions, f[2], alts)),
    }
}

/// The rows of a default-format output (`iter_rows`).
pub struct Rows<'a> {
    path: String,
    dataset: String,
    keys: RowKeys,
    specs: Vec<FamilySpec>,
    substitutions: Option<&'a Substitutions>,
    lines: Option<TextLines>,
    failure: Option<Failure>,
}

/// Every row of at least seven columns of a default-format reference, a line that starts
/// with `#` being none; a file that cannot be opened or read is refused with its path, as
/// `vcf_lines` refuses a VCF. `family_specs` are the specs of `--family-columns`, none for
/// rows without family classes. With `substitutions` (the reference's input, as
/// `substitution_alts` reads it) a row whose input record holds its allele as a
/// substitution is marked so (`holds_allele`), which its variant class is read with.
pub fn iter_rows<'a, S: AsRef<str>>(
    path: &str,
    dataset: &str,
    family_specs: &[S],
    substitutions: Option<&'a Substitutions>,
) -> Rows<'a> {
    let (lines, failure) = match opener(path) {
        Ok(lines) => (Some(lines), None),
        Err(error) => (None, Some(not_readable_text(path, &error))),
    };
    Rows {
        path: path.to_string(),
        dataset: dataset.to_string(),
        keys: RowKeys {
            strand: Finder::new("STRAND"),
            flags: Finder::new("FLAGS"),
        },
        specs: family_specs
            .iter()
            .map(|spec| FamilySpec::read(spec.as_ref()))
            .collect(),
        substitutions: substitutions.filter(|held| !held.is_empty()),
        lines,
        failure,
    }
}

impl Rows<'_> {
    /// The next row as the reader's line holds it. The iterator gives the same rows with
    /// members of their own.
    pub fn next_view(&mut self) -> Option<Result<RowView<'_>, Failure>> {
        if let Some(failure) = self.failure.take() {
            return Some(Err(failure));
        }
        loop {
            let lines = self.lines.as_mut()?;
            match lines.advance() {
                Ok(true) => {}
                Ok(false) => {
                    self.lines = None;
                    return None;
                }
                Err(error) => {
                    self.lines = None;
                    return Some(Err(not_readable_text(&self.path, &error)));
                }
            }
            let line = lines.line();
            let tabs = memchr::memchr_iter(b'\t', line.as_bytes()).take(6).count();
            if tabs == 6 && !line.starts_with('#') {
                break;
            }
        }
        let line = self.lines.as_ref()?.line();
        Some(Ok(row_view(
            line,
            &self.dataset,
            &self.keys,
            &self.specs,
            self.substitutions,
        )))
    }
}

impl Iterator for Rows<'_> {
    type Item = Result<Row, Failure>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.next_view()?.map(|view| view.to_row()))
    }
}

/// The rows of a default-format output, one input record at a time (`iter_records`).
pub struct Records<'a> {
    rows: Rows<'a>,
    first_of_next: Option<Row>,
}

/// The rows of one input record at a time: the consecutive rows of a reference sharing a
/// `record_key`, which Ensembl VEP writes together (every row of a variant comes out of one
/// output-factory call); `substitutions` is `iter_rows`'s. A file that cannot be read on is
/// refused in place of the record the reading stopped in.
pub fn iter_records<'a>(
    path: &str,
    dataset: &str,
    substitutions: Option<&'a Substitutions>,
) -> Records<'a> {
    Records {
        rows: iter_rows::<&str>(path, dataset, &[], substitutions),
        first_of_next: None,
    }
}

impl Iterator for Records<'_> {
    type Item = Result<Vec<Row>, Failure>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut group: Vec<Row> = self.first_of_next.take().into_iter().collect();
        loop {
            let view = match self.rows.next_view() {
                Some(Ok(view)) => view,
                Some(Err(failure)) => return Some(Err(failure)),
                None if group.is_empty() => return None,
                None => return Some(Ok(group)),
            };
            if group
                .first()
                .is_some_and(|first| first.record_key() != view.record_key())
            {
                self.first_of_next = Some(view.to_row());
                return Some(Ok(group));
            }
            group.push(view.to_row());
        }
    }
}

/// The rows of one input record, held by every sample that keeps the record.
pub type RecordRows = Arc<[Row]>;

/// One step of a bounded sample (`reservoir_step` with its generator): the bucket, the
/// item offered, the count of items offered to the bucket so far, the bucket's capacity.
pub type ReservoirStep<'s> = dyn FnMut(&mut Vec<RecordRows>, RecordRows, usize, usize) + 's;

/// What one pass over the references yields for a record-class key.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RecordClassPass {
    /// The classes of every classed record, by (dataset id, Uploaded_variation, Location).
    pub classes_of: HashMap<(String, String, String), BTreeSet<&'static str>>,
    /// Per class, a bounded uniform sample of the records that show it, each as its rows.
    pub samples: BTreeMap<&'static str, Vec<RecordRows>>,
    /// The count of records per class.
    pub seen: BTreeMap<&'static str, usize>,
}

/// One streaming pass over the references for the record classes of `key`, each record
/// classed from its rows and its input record's alleles. A classed record is offered to the
/// sample of each of its classes in the classes' sorted order through `step`, which holds
/// the generator `random.Random(seed)` of the pass, with `RECORD_CLASS_RESERVOIR` as the
/// capacity.
pub fn record_class_pass(
    datasets: &[Dataset],
    key: &str,
    step: &mut ReservoirStep<'_>,
) -> Result<RecordClassPass, Failure> {
    let Some(classes_under_key) = record_class_function(key) else {
        return Err(Failure::own(format!(
            "{} is no record-class key",
            py_repr_str(key)
        )));
    };
    let mut pass = RecordClassPass::default();
    for (sid, reference, input) in datasets {
        let alleles = input_alleles(input)?;
        let substitutions = substitution_alts(input)?;
        for rows in iter_records(reference, sid, Some(&substitutions)) {
            let rows = rows?;
            let (name, location) = rows[0].record_key();
            let chromosome = location.split_once(':').map_or(location, |(head, _)| head);
            let pairs = alleles.get(chromosome, name).unwrap_or_default();
            let classes = classes_under_key(&rows, pairs);
            if classes.is_empty() {
                continue;
            }
            let record = (sid.clone(), name.to_string(), location.to_string());
            let rows: RecordRows = rows.into();
            for class in &classes {
                let seen = pass.seen.entry(*class).or_default();
                *seen += 1;
                let bucket = pass.samples.entry(*class).or_default();
                step(bucket, Arc::clone(&rows), *seen, RECORD_CLASS_RESERVOIR);
            }
            pass.classes_of.insert(record, classes);
        }
    }
    Ok(pass)
}

#[cfg(test)]
#[path = "rows_tests.rs"]
mod tests;
