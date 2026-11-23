// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The input index and the VCF parser.

use crate::entry::{Entry, Fields};
use crate::error::{Error, Result};
use crate::input::{dedup_keys, line_slices};
use crate::pyfmt;
use crate::rules::record_alleles::vcf_alleles_as_spelt;
use crate::rules::structural::structural_line_location;
use crate::self_parity::reverse_complement;
use crate::shared::{
    allele_is_indel, integer_of, is_structural_allele, reference_allele_is_empty,
    replaced_reference, sequence_location,
};
use indexmap::IndexMap;
use std::collections::VecDeque;

/// The strand an Ensembl-format input line's fifth field names.
pub const ENSEMBL_FORMAT_STRANDS: &[(&str, i64)] =
    &[("+", 1), ("1", 1), ("+1", 1), ("-", -1), ("-1", -1)];

/// The fields a record entry takes from a line's first seven columns, in column order.
const FIXED_COLUMNS: [&str; 7] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"];

/// Input records by identity, for mapping an output record back to its ordinal.
/// Identical input lines are matched in order of appearance.
#[derive(Debug)]
pub struct InputIndex {
    pub unmatched: IndexMap<String, VecDeque<i64>>,
}

impl InputIndex {
    pub fn new(inputs: &[String], key: &mut dyn FnMut(&str) -> String) -> Self {
        let mut unmatched: IndexMap<String, VecDeque<i64>> = IndexMap::new();
        for (ordinal, line) in (0i64..).zip(inputs) {
            unmatched.entry(key(line)).or_default().push_back(ordinal);
        }
        InputIndex { unmatched }
    }

    pub fn record(&mut self, key: &str) -> Option<i64> {
        self.unmatched.get_mut(key)?.pop_front()
    }
}

/// A VCF data line's identity: CHROM, POS, ID, REF and ALT.
pub fn vcf_line_key(line: &str) -> String {
    line_key(line).to_string()
}

/// `vcf_line_key` as a part of `line`: the text before its fifth tab, the whole line where
/// it has fewer.
pub(crate) fn line_key(line: &str) -> &str {
    match memchr::memchr_iter(b'\t', line.as_bytes()).nth(4) {
        Some(end) => &line[..end],
        None => line,
    }
}

/// What the comparison reads of an input line, a part of `line`: the text before its fifth
/// tab, which holds its key (`line_key`) and the five columns `input_line_record` and
/// `RecordRows` read, lengthened to the end of its fifth field between white space where
/// that field ends later (`ensembl_format_record` reads five such fields, and an empty or
/// blank column among the first five puts the fifth in a later column). A line of fewer
/// than five tabs or fewer than five fields is whole. A reader that wants more of a line
/// than these lengthens this part: `input_records` keeps no more.
pub(crate) fn record_prefix(line: &str) -> &str {
    let Some(fifth_tab) = memchr::memchr_iter(b'\t', line.as_bytes()).nth(4) else {
        return line;
    };
    let (mut fields, mut inside) = (0, false);
    for (at, c) in line.char_indices() {
        if !pyfmt::char_is_space(c) {
            fields += usize::from(!inside);
            inside = true;
        } else if inside {
            inside = false;
            if fields == 5 {
                return &line[..at.max(fifth_tab)];
            }
        }
    }
    line
}

/// `int(digits)` of an input line's position that `str.isdigit` has passed. `int` reads no
/// digit outside the decimal ones (a superscript) and no number longer than the interpreter
/// converts: the `ValueError` ends the run. A number past `i64` ends the run in the same
/// words.
fn position_of(digits: &str) -> i64 {
    if let Some(position) = pyfmt::int_of(digits) {
        return position;
    }
    match pyfmt::int_checked(digits) {
        Ok(_) => crate::uncaught!("{}", pyfmt::IntError::Invalid(pyfmt::repr_str(digits))),
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// The chromosome, first position, reference allele and alternate alleles of an input line
/// in the Ensembl default format: the chromosome, an integer start and end, the allele string
/// `REF/ALT[/ALT...]` with `-` for an empty allele, and the strand, separated by white
/// space. The alleles come in upper case on the forward strand, those of a line on the reverse
/// strand reverse-complemented. An insertion written alone has its start one above its end and
/// an empty reference allele. `None` for a line of any other shape.
pub fn ensembl_format_record(line: &str) -> Option<(String, i64, String, Vec<String>)> {
    let fields = pyfmt::split_whitespace(line);
    if fields.len() < 5 || !(pyfmt::is_digit(fields[1]) && pyfmt::is_digit(fields[2])) {
        return None;
    }
    if !fields[3].contains('/') {
        return None;
    }
    let (_, strand) = ENSEMBL_FORMAT_STRANDS
        .iter()
        .find(|(spelt, _)| *spelt == fields[4])?;
    let mut alleles = fields[3].split('/').map(|allele| {
        let allele = if allele == "-" {
            String::new()
        } else {
            pyfmt::upper(allele)
        };
        if *strand < 0 {
            reverse_complement(&allele)
        } else {
            allele
        }
    });
    let reference = alleles.next().unwrap_or_default();
    let alternates: Vec<String> = alleles.collect();
    Some((
        fields[0].to_string(),
        position_of(fields[1]),
        reference,
        alternates,
    ))
}

/// The chromosome, first position, reference allele and alternate alleles an input line
/// gives, the alleles in upper case on the forward strand. The line's fourth field tells its
/// format: a VCF line's REF is letters alone, beside an integer POS and an ALT column, and an
/// Ensembl-format line's allele string holds a `/` (`ensembl_format_record`). `None`
/// for a line of neither shape.
pub fn input_line_record(line: &str) -> Option<(String, i64, String, Vec<String>)> {
    let columns: Vec<&str> = line.splitn(6, '\t').collect();
    if columns.len() >= 5
        && pyfmt::is_digit(columns[1])
        && !columns[3].is_empty()
        && columns[3]
            .bytes()
            .all(|letter| letter.is_ascii_alphabetic())
    {
        return Some((
            columns[0].to_string(),
            position_of(columns[1]),
            pyfmt::upper(columns[3]),
            pyfmt::upper(columns[4])
                .split(',')
                .map(str::to_string)
                .collect(),
        ));
    }
    ensembl_format_record(line)
}

/// The ordinal of the Ensembl-format input line a VCF output line is written from, for a
/// line whose own columns name no input record: an engine writes one VCF line per input line,
/// in order, at the input line's start or, where it pads the alleles with the base before
/// them, one position below. `ordinal` is the output line's own, `chromosome` and
/// `position` its CHROM and POS; `None` where the input line there is no such line.
pub fn ensembl_line_ordinal(
    inputs: &[String],
    ordinal: i64,
    chromosome: &str,
    position: &str,
) -> Option<i64> {
    let line = usize::try_from(ordinal)
        .ok()
        .and_then(|at| inputs.get(at))?;
    let (named, start, _, _) = ensembl_format_record(line)?;
    if named != chromosome || !pyfmt::is_digit(position) {
        return None;
    }
    let position = match pyfmt::int_of(position) {
        Some(position) => position,
        None => match pyfmt::int_checked(position) {
            // A number past `i64` is neither the line's start nor the position below it.
            Ok(_) => return None,
            Err(error) => crate::uncaught!("{error}"),
        },
    };
    (position == start || position == start - 1).then_some(ordinal)
}

/// The CSQ field names of a `##INFO=<ID=CSQ` header line's `Format:` list.
pub fn csq_format_fields(line: &str) -> Result<Vec<String>> {
    let Some((_, mut format)) = line.split_once("Format: ") else {
        return Err(Error::comparator(format!(
            "CSQ header without a Format list: {line}"
        )));
    };
    while let Some(shorter) = format.strip_suffix("\">") {
        format = shorter;
    }
    Ok(format.split('|').map(str::to_string).collect())
}

/// The columns of a default or tab output's `#` header line, empty without one.
pub fn header_columns(text: &str) -> Vec<String> {
    for line in line_slices(text) {
        if line.starts_with('#') && !line.starts_with("##") {
            return line[1..].split('\t').map(str::to_string).collect();
        }
    }
    Vec::new()
}

/// The POS of a VCF line and the Location the default format writes for its rows: the
/// line's span, without the first base where Parser/VCF.pm:300-336 drop it
/// (`vcf_alleles_as_spelt`), an emptied reference allele's two neighbouring positions.
/// `None` for a POS that is not a run of digits.
///
/// A span that ends past `i64` is refused as a POS past `i64` is (`integer_of`).
pub fn vcf_line_span(
    chromosome: &str,
    position: &str,
    reference: &str,
    alternates: &[String],
) -> Result<Option<(i64, String)>> {
    if position.is_empty() || !position.bytes().all(|digit| digit.is_ascii_digit()) {
        return Ok(None);
    }
    let first = integer_of(position, "VCF line POS")?;
    let dropped = vcf_alleles_as_spelt(reference, alternates).0 != reference;
    let bound = |offset: i128| {
        let bound = i128::from(first) + offset;
        i64::try_from(bound).or_else(|_| integer_of(&bound.to_string(), "VCF line POS"))
    };
    let start = bound(i128::from(dropped))?;
    let end = bound(pyfmt::len(reference) as i128 - 1)?;
    Ok(Some((first, sequence_location(chromosome, start, end))))
}

/// The VCF format: one entry per CSQ entry, keyed by the input record ordinal, the
/// allele, the feature and the feature type, plus one record entry per line carrying the
/// seven fixed columns, the INFO without CSQ and whether a CSQ was present. A line
/// matching no input record is keyed by its output ordinal; its CSQ entries carry the ordinal
/// of the Ensembl-format input line it is written from, where it is (`Entry.line`). A CSQ
/// entry carries the reference bases its allele stands in for and the Location of the
/// default format (`Entry.reference`, `Entry.allele_location`) where the line's POS, REF
/// and ALT give the allele at the line's span. A structural allele's entry stands in for
/// no bases and carries the Location of its line on the structural path: the base after
/// POS to the end INFO/SVLEN or INFO/END gives (`structural_line_location`).
pub fn parse_vcf(text: &str, inputs: &[String]) -> Result<Vec<Entry>> {
    vcf_entries(line_slices(text), None, inputs)
}

/// `parse_vcf` of some lines of an output, as `text_lines` gives them: its `#` lines and
/// the data lines a comparison reads, in file order. `ordinals` holds one number a data
/// line of `lines`, in order: the line's ordinal among the data lines of the whole output,
/// which labels a line that matches no input record. A `#` line has no ordinal.
pub fn parse_vcf_lines(
    lines: &[String],
    ordinals: &[i64],
    inputs: &[String],
) -> Result<Vec<Entry>> {
    vcf_entries(lines.iter().map(String::as_str), Some(ordinals), inputs)
}

fn vcf_entries<'a>(
    lines: impl Iterator<Item = &'a str>,
    ordinals: Option<&[i64]>,
    inputs: &[String],
) -> Result<Vec<Entry>> {
    let mut csq_fields: Vec<String> = Vec::new();
    let mut entries: Vec<Entry> = Vec::new();
    let mut index = InputIndex::new(inputs, &mut vcf_line_key);
    // A whole output numbers its own data lines.
    let mut counted = 0i64..;
    let mut given = ordinals.map(|ordinals| ordinals.iter().copied());
    for line in lines {
        if line.starts_with("##INFO=<ID=CSQ") {
            csq_fields = csq_format_fields(line)?;
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let mut columns = line.splitn(9, '\t');
        let mut cols = [""; 8];
        for column in &mut cols {
            *column = columns.next().ok_or_else(|| {
                Error::comparator(format!("VCF line has fewer than 8 columns: {line}"))
            })?;
        }
        let record = index.record(line_key(line));
        let ordinal = match &mut given {
            None => counted.next(),
            Some(given) => given.next(),
        }
        .ok_or_else(|| Error::comparator("a data line has no ordinal in its output"))?;
        let rec_label = match record {
            Some(record) => format!("rec{record}"),
            None => format!("line{ordinal}"),
        };
        let input_line = match record {
            Some(_) => record,
            None => ensembl_line_ordinal(inputs, ordinal, cols[0], cols[1]),
        };
        let mut info_rest: Vec<&str> = Vec::new();
        let mut csq: Option<&str> = None;
        for item in cols[7].split(';') {
            match item.strip_prefix("CSQ=") {
                Some(value) => csq = Some(value),
                None => info_rest.push(item),
            }
        }
        let mut fields = Fields::with_capacity(FIXED_COLUMNS.len() + 2);
        for (name, value) in FIXED_COLUMNS.iter().zip(cols) {
            fields.insert(name.to_string(), value.to_string());
        }
        fields.insert("INFO_without_CSQ".to_string(), info_rest.join(";"));
        let csq_present = if csq.is_some() { "yes" } else { "no" };
        fields.insert("CSQ_present".to_string(), csq_present.to_string());
        entries.push(Entry {
            key: format!("{rec_label}|record"),
            record,
            fields,
            ..Entry::default()
        });
        let Some(csq) = csq else {
            continue;
        };
        let reference_allele = cols[3];
        let alternate_alleles: Vec<String> = cols[4].split(',').map(str::to_string).collect();
        let alternates: Vec<&str> = cols[4].split(',').collect();
        let span = vcf_line_span(cols[0], cols[1], reference_allele, &alternate_alleles)?;
        // Read at the line's first structural allele: most lines carry none.
        let mut structural_span: Option<Option<String>> = None;
        for chunk in csq.split(',') {
            let width = memchr::memchr_iter(b'|', chunk.as_bytes()).count() + 1;
            if width != csq_fields.len() {
                return Err(Error::comparator(format!(
                    "CSQ chunk width on {rec_label}: {chunk}"
                )));
            }
            let mut cfields = Fields::with_capacity(width);
            for (name, value) in csq_fields.iter().zip(chunk.split('|')) {
                cfields.insert(name.clone(), value.to_string());
            }
            let field = |name: &str| cfields.get(name).cloned().unwrap_or_default();
            let allele = field("Allele");
            let feature = field("Feature");
            let reference = match &span {
                Some((first, location)) => replaced_reference(
                    *first,
                    reference_allele,
                    &alternates,
                    Some(location),
                    &allele,
                )?,
                None => None,
            };
            let allele_location = match (&span, &reference) {
                (Some((_, location)), Some(_)) => Some(location.clone()),
                (_, None) if is_structural_allele(&allele) => structural_span
                    .get_or_insert_with(|| {
                        structural_line_location(
                            cols[0],
                            cols[1],
                            reference_allele,
                            cols[4],
                            &info_rest.join(";"),
                        )
                    })
                    .clone(),
                _ => None,
            };
            entries.push(Entry {
                key: format!("{rec_label}|{allele}|{feature}|{}", field("Feature_type")),
                record,
                location: None,
                consequence: field("Consequence"),
                indel: allele_is_indel(reference_allele, &alternates, &allele, false),
                insertion: Some(reference_allele_is_empty(
                    reference_allele,
                    &alternates,
                    false,
                )),
                reference,
                line: input_line,
                allele_location,
                allele,
                feature,
                fields: cfields,
            });
        }
    }
    dedup_keys(&mut entries);
    Ok(entries)
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_part_kept_of_an_input_line_reads_as_the_line_reads() {
        for (line, kept) in [
            // The columns after the fifth are not kept.
            (
                "21\t5\trs1\tA\tT\t.\tPASS\tAC=1;AF=0.5 x",
                "21\t5\trs1\tA\tT",
            ),
            ("21\t5\trs1\tA\tT\t", "21\t5\trs1\tA\tT"),
            // A line of fewer than five tabs is whole.
            ("21\t5\trs1\tA\tT", "21\t5\trs1\tA\tT"),
            ("21 5 5 A/T + rs1 more", "21 5 5 A/T + rs1 more"),
            ("one", "one"),
            ("", ""),
            // An empty or blank column among the first five puts the fifth field between
            // white space in a later column, which is kept to its end.
            ("21\t100\t\t101\tA/T\t+\t99\tAC=1", "21\t100\t\t101\tA/T\t+"),
            (
                "21\t100\t \t101\tA/T\t-1 9\tAC=1",
                "21\t100\t \t101\tA/T\t-1",
            ),
            ("\t21\t100\t101\tA/T\t\t+\tAC=1", "\t21\t100\t101\tA/T\t\t+"),
            // A field that ends the line, and a line of fewer than five fields, are whole.
            ("21\t100\t\t101\tA/T\t+", "21\t100\t\t101\tA/T\t+"),
            ("21\t\t\t\t\tx", "21\t\t\t\t\tx"),
            ("\t\t\t\t\t", "\t\t\t\t\t"),
            // White space of any script parts the fields.
            (
                "21\u{a0}100\t101\tA/T\t+\tx\ty z",
                "21\u{a0}100\t101\tA/T\t+\tx",
            ),
        ] {
            let part = record_prefix(line);
            assert_eq!(part, kept, "{line:?}");
            assert_eq!(vcf_line_key(part), vcf_line_key(line), "{line:?}");
            assert_eq!(input_line_record(part), input_line_record(line), "{line:?}");
            assert_eq!(
                ensembl_format_record(part),
                ensembl_format_record(line),
                "{line:?}"
            );
            assert_eq!(
                part.splitn(6, '\t').take(5).collect::<Vec<_>>(),
                line.splitn(6, '\t').take(5).collect::<Vec<_>>(),
                "{line:?}"
            );
        }
        // The line of the third shape is an Ensembl-format record only by its sixth column.
        assert_eq!(
            ensembl_format_record("21\t100\t\t101\tA/T\t+\t99"),
            Some((
                "21".to_string(),
                100,
                "A".to_string(),
                vec!["T".to_string()]
            ))
        );
    }
}
