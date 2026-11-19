// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` source over a tabix-indexed VCF, read as
//! `AnnotationSource/File/VCF.pm` reads one through ensembl-io's
//! `Parser/BaseVCF4.pm`.
//!
//! A record's coordinates follow ensembl-io: its start moves one past `POS`
//! for a structural record and for an indel, its end comes from `SVLEN`,
//! `END`, the structural start or the reference's length (`BaseVCF4.pm:255-317`).
//! A structural record, or any record against a structural variant, matches
//! through the base predicates ([`base_overlap`]); a sequence record matches an
//! `overlap` source over `POS..POS + length(REF) - 1` and an `exact` source
//! allele by allele, each pair of reference and alternate trimmed to its
//! minimal form before the comparison (ensembl-variation
//! `Utils/Sequence.pm:1144-1245`); `within` and `surrounding` never match a
//! sequence record (`VCF.pm:389-422`). The INFO fields the source names,
//! the record's FILTER entries, its first ID or its coordinates, and its
//! QUAL become the output records (`VCF.pm:173-318`).

use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use bstr::ByteSlice;
use indexmap::IndexMap;
use noodles::bgzf;
use vep_core::custom_annotation::{CustomFieldValue, CustomRecord, CustomValue};
use vep_core::sv_routing::{sv_so_term_and_class, sv_type_abbreviation};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{first_vcf_id, InputVariant};
use vep_core::vcf_header::fileformat_refusal;

use super::grammar::{perl_number, perl_split, truthy, CustomSpec, MatchType};
use super::source::{
    base_overlap, tabix_reason, whole_number, Coordinate, CustomError, CustomSource, FormatRules,
    RecordGroup, TabixSource, Verdict,
};
use super::CustomRunOptions;
use crate::tabix::allele_match::{dashed, trim_directions, trim_sequences};
use crate::tabix::{lossy_text, BoundedLines, TabixRecord};
use crate::traits::PrefetchData;

/// The `##` lines of a VCF source, as ensembl-io reads them at open
/// (`BaseVCF4.pm:86-163`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct VcfMetadata {
    /// The value of the last `##source=` line.
    pub source: Option<String>,
    /// Each `##INFO` line's `Description` by its `ID`, where the line gives
    /// both; a later line with the same `ID` replaces the earlier.
    pub info_descriptions: IndexMap<String, String>,
}

impl VcfMetadata {
    /// Reads the `#` lines that open the bgzipped file at `path`, each less its
    /// `\n` or `\r\n` and nothing more, as htslib hands the header lines to
    /// ensembl-io (`VCF4Tabix.pm:55-58`), so a second `\r` is part of the
    /// line's value. A `##fileformat` line naming a VCF major version other
    /// than 4, or a value with no `<digits>.<digits>` version, refuses the file
    /// with Ensembl VEP's text behind the line read, as ensembl-io's parser
    /// ends the run on it while the source opens (`BaseVCF4.pm:125, 129`).
    pub fn read(path: &Path) -> Result<Self, CustomError> {
        let open_error = |reason: String| CustomError::Open {
            file: path.display().to_string(),
            reason,
        };
        let file = std::fs::File::open(path).map_err(|e| open_error(e.to_string()))?;
        let reader = bgzf::io::Reader::new(BufReader::new(file));
        let mut metadata = VcfMetadata::default();
        for line in BoundedLines::new(reader, path) {
            let line = line.map_err(|e| open_error(tabix_reason(e)))?;
            if !line.starts_with(b"#") {
                break;
            }
            metadata.read_line(&lossy_text(line)).map_err(open_error)?;
        }
        Ok(metadata)
    }

    /// The metadata of the given header lines, or the refusal of a
    /// `##fileformat` line among them.
    pub fn parse_lines<'a>(lines: impl IntoIterator<Item = &'a str>) -> Result<Self, String> {
        let mut metadata = VcfMetadata::default();
        for line in lines {
            metadata.read_line(line)?;
        }
        Ok(metadata)
    }

    /// `read_metadata` (`BaseVCF4.pm:86-100`): a line `##<key>=<value>`,
    /// blanks allowed after the `##`, is a key of word characters and a
    /// non-empty value. The error is the refusal of a `##fileformat` line
    /// ([`fileformat_refusal`]), the line quoted before Ensembl VEP's text.
    fn read_line(&mut self, line: &str) -> Result<(), String> {
        if let Some(refusal) = fileformat_refusal(line.as_bytes()) {
            return Err(format!("{line}: {refusal}"));
        }
        let Some(rest) = line.strip_prefix("##") else {
            return Ok(());
        };
        let rest = rest.trim_start_matches(|c: char| c.is_ascii_whitespace() || c == '\u{0B}');
        let key_end = rest
            .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
            .unwrap_or(rest.len());
        let Some(value) = rest[key_end..].strip_prefix('=') else {
            return Ok(());
        };
        if key_end == 0 || value.is_empty() {
            return Ok(());
        }
        match &rest[..key_end] {
            "INFO" => {
                let entry = structured_entry(value);
                if let (Some(Some(id)), Some(Some(description))) =
                    (entry.get("ID"), entry.get("Description"))
                {
                    self.info_descriptions
                        .insert(id.clone(), description.clone());
                }
            }
            "source" => self.source = Some(value.to_string()),
            _ => {}
        }
        Ok(())
    }
}

/// `_parse_metadata_line` on a structured line (`BaseVCF4.pm:133-152`): the
/// angle brackets go, the commas inside the last quoted text are kept out of
/// the split, each `,`-separated item is a key and the piece after its first
/// `=` (absent for a key alone), and a value loses every `"`.
fn structured_entry(value: &str) -> IndexMap<String, Option<String>> {
    let mut data: String = value.chars().filter(|c| !matches!(c, '<' | '>')).collect();
    let quotes: Vec<usize> = data.match_indices('"').map(|(i, _)| i).collect();
    if let Some(&last) = quotes.last() {
        // `(.*)(".+")(.*)`: the group opens at the last quote that has a
        // quote at least two characters after it and closes at the last quote.
        if let Some(&first) = quotes.iter().rev().find(|&&q| q + 2 <= last) {
            let protected = data[first..=last].replace(',', "!#!");
            data.replace_range(first..=last, &protected);
        }
    }
    let mut entry = IndexMap::new();
    for item in perl_split(&data, ',') {
        // `my ($key, $value) = split('=', $meta)`: the value is the text between
        // the first and second `=`, an empty one kept, none without an `=`.
        let (key, value) = match item.split_once('=') {
            Some((key, rest)) => (key, Some(rest.split('=').next().unwrap_or(""))),
            None => (item, None),
        };
        let value = value.map(|v| v.replace('"', "").replace("!#!", ","));
        entry.insert(key.to_string(), value);
    }
    entry
}

/// Whether the output reads the file's `##INFO` descriptions for its header:
/// the name ends in `.vcf`, `.vcf.gz` or `.vcf.bgz` in any case
/// (`OutputFactory/BaseTab.pm:86`).
fn vcf_named(file: &str) -> bool {
    let lower = file.to_ascii_lowercase();
    lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz")
}

/// A `--custom` source of `format=vcf`.
pub struct VcfSource {
    base: TabixSource,
    options: CustomRunOptions,
    metadata: VcfMetadata,
    descriptions_apply: bool,
    /// The INFO keys read from every candidate record.
    keys: InfoKeys,
}

impl VcfSource {
    pub fn open(
        spec: CustomSpec,
        synonyms: Arc<ChromosomeSynonyms>,
        options: CustomRunOptions,
    ) -> Result<Self, CustomError> {
        let base = TabixSource::open(spec, synonyms)?;
        let metadata = VcfMetadata::read(Path::new(&base.spec().file))?;
        let descriptions_apply = vcf_named(&base.spec().file);
        let keys = InfoKeys::new(base.spec().fields.as_deref().unwrap_or(&[]));
        Ok(VcfSource {
            base,
            options,
            metadata,
            descriptions_apply,
            keys,
        })
    }

    /// The file's `##` metadata.
    pub fn metadata(&self) -> &VcfMetadata {
        &self.metadata
    }

    /// `_get_record_name` (`VCF.pm:283-297`): the first `;`-separated ID
    /// ([`first_vcf_id`]), or `<chromosome>:<start>-<end>` under `coords`,
    /// when the ID column is empty, Perl's false or `.`.
    fn record_name(&self, candidate: &Candidate<'_>) -> String {
        match first_vcf_id(candidate.line.id) {
            Some(id) if !self.base.spec().report_coords && truthy(id) && id != "." => {
                id.to_string()
            }
            _ => format!(
                "{}:{}-{}",
                candidate.line.chr,
                candidate.start.whole(),
                candidate.end.whole()
            ),
        }
    }

    /// The data of the fields the source names, for one record (`VCF.pm:181-224`).
    fn fields_data(&self, candidate: &Candidate<'_>) -> Vec<(String, FieldData)> {
        let Some(fields) = &self.base.spec().fields else {
            return Vec::new();
        };
        let alt_count = candidate.alts.len();
        let clinvar = self.metadata.source.as_deref() == Some("ClinVar");
        let mut data: Vec<(String, FieldData)> = Vec::new();
        for (field, value) in fields.iter().zip(&candidate.info.fields) {
            let Some(value) = value else {
                continue;
            };
            // A key without a value is a flag, written as `1` (`:193`).
            let value = value.unwrap_or("1").to_string();
            let data_value = if !(self.options.multi_allelic || clinvar) && value.contains(',') {
                let mut split = perl_split(&value, ',');
                // A leading value for the reference is dropped (`:203`).
                if split.len() == alt_count + 1 {
                    split.remove(0);
                }
                if split.len() == alt_count {
                    FieldData::PerAlt(split.into_iter().map(str::to_string).collect())
                } else {
                    FieldData::Whole(value)
                }
            } else {
                FieldData::Whole(value)
            };
            data.push((field.clone(), data_value));
        }
        let names_filter = fields.iter().any(|f| f == "FILTER");
        if names_filter || !self.options.suppress_filter {
            // `.=` onto a named INFO key of that name (`:221`).
            let mut filter = match data.iter().position(|(name, _)| name == "FILTER") {
                Some(index) => match data.remove(index).1 {
                    FieldData::Whole(text) => text,
                    FieldData::PerAlt(_) | FieldData::Filter(_) => String::new(),
                },
                None => String::new(),
            };
            filter.push_str(candidate.line.filter);
            data.push((
                "FILTER".to_string(),
                FieldData::Filter(
                    perl_split(&filter, ';')
                        .into_iter()
                        .map(str::to_string)
                        .collect(),
                ),
            ));
        }
        data
    }

    /// `_create_records` (`VCF.pm:173-266`): one record per matched allele
    /// under `exact`, carrying the allele and each per-allele value at the
    /// matched alternate's index; else one record for the whole variant, a
    /// per-allele value rejoined with `,` in alternate order.
    fn create_records(
        &self,
        candidate: &Candidate<'_>,
        matches: Option<&[(usize, &InputVariant)]>,
        scores: bool,
    ) -> Vec<CustomRecord> {
        let fields_data = self.fields_data(candidate);
        let name = self.record_name(candidate);
        let score = if scores { candidate.line.score() } else { None };
        let fields = |index: Option<usize>| -> IndexMap<String, CustomFieldValue> {
            fields_data
                .iter()
                .map(|(field, data)| {
                    let value = match data {
                        FieldData::Whole(text) => CustomFieldValue::Text(text.clone()),
                        FieldData::Filter(entries) => CustomFieldValue::List(entries.clone()),
                        FieldData::PerAlt(values) => CustomFieldValue::Text(match index {
                            Some(index) => values[index].clone(),
                            None => values.join(","),
                        }),
                    };
                    (field.clone(), value)
                })
                .collect()
        };
        match matches {
            Some(matches) => matches
                .iter()
                .map(|(index, variant)| CustomRecord {
                    name: name.clone(),
                    allele: Some(variant.alt_allele().to_vec()),
                    fields: fields(Some(*index)),
                    score: score.clone(),
                })
                .collect(),
            None => vec![CustomRecord {
                name,
                allele: None,
                fields: fields(None),
                score,
            }],
        }
    }
}

/// One record of the file as a candidate for an input record: its columns,
/// the INFO values it is read for, its alternates and its feature
/// coordinates.
struct Candidate<'a> {
    line: VcfLine<'a>,
    info: InfoValues<'a>,
    alts: Vec<&'a str>,
    /// `get_start`.
    start: Coordinate,
    /// `get_end`.
    end: Coordinate,
}

impl<'a> Candidate<'a> {
    /// `pos` is the record's `POS` and `end_of` reads its `END`: each as Perl
    /// reads the cell ([`TabixSource::start_of`], [`perl_end`]) or, for an
    /// `exact` source, in whole numbers ([`TabixSource::keyed_start_of`],
    /// [`whole_end`]).
    fn of(
        record: &'a TabixRecord,
        keys: &InfoKeys,
        pos: Coordinate,
        end_of: fn(&str) -> Coordinate,
    ) -> Self {
        let line = VcfLine {
            pos,
            ..VcfLine::of(record)
        };
        let info = line.info_values(keys);
        let alts = line.alternatives();
        let structural = line.structural_shape();
        let start = line.feature_start(&alts, structural);
        let end = line.feature_end(start, &info, structural, end_of);
        Candidate {
            line,
            info,
            alts,
            start,
            end,
        }
    }
}

/// The INFO keys a source reads from every candidate record, each in the
/// delimited form `;KEY` an item takes inside the column, built once at
/// open. A key holding `;` or `=`, the column's own separators, names no
/// item: its needle is empty.
struct InfoKeys {
    svlen: Vec<u8>,
    end: Vec<u8>,
    svtype: Vec<u8>,
    /// One needle per field the source names, in the source's order.
    fields: Vec<Vec<u8>>,
}

impl InfoKeys {
    fn new(fields: &[String]) -> Self {
        let needle = |key: &str| -> Vec<u8> {
            if key.contains([';', '=']) {
                return Vec::new();
            }
            let mut needle = Vec::with_capacity(key.len() + 1);
            needle.push(b';');
            needle.extend_from_slice(key.as_bytes());
            needle
        };
        InfoKeys {
            svlen: needle("SVLEN"),
            end: needle("END"),
            svtype: needle("SVTYPE"),
            fields: fields.iter().map(|field| needle(field)).collect(),
        }
    }
}

/// The INFO values a candidate is read for, each as `get_info` holds the key
/// (`BaseVCF4.pm:576-589`): `None` for a key the column lacks, `Some(None)`
/// for a key alone, else the text between the key's first and second `=`,
/// empty for `KEY=`; a later item with the key replaces an earlier. The keys
/// asked for are read from the column by [`VcfLine::info_values`], in place
/// of the hash Ensembl VEP builds over every key.
struct InfoValues<'a> {
    svlen: Option<Option<&'a str>>,
    end: Option<Option<&'a str>>,
    svtype: Option<Option<&'a str>>,
    /// One entry per requested field, in the source's order.
    fields: Vec<Option<Option<&'a str>>>,
}

/// The data of one named field before the records are built
/// (`$fields_data`, `VCF.pm:181-224`).
enum FieldData {
    /// One value for every allele.
    Whole(String),
    /// One value per alternate allele, in the ALT column's order.
    PerAlt(Vec<String>),
    /// The FILTER column's `;`-separated entries.
    Filter(Vec<String>),
}

impl CustomSource for VcfSource {
    fn spec(&self) -> &CustomSpec {
        self.base.spec()
    }

    fn chromosomes(&self) -> Vec<String> {
        self.base.chromosomes()
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
        self.base.prefetch(variants)
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), CustomError> {
        self.base.annotate(variants, data, self)
    }

    fn field_descriptions(&self) -> Option<&IndexMap<String, String>> {
        self.descriptions_apply
            .then_some(&self.metadata.info_descriptions)
    }
}

impl FormatRules for VcfSource {
    /// `_record_overlaps_VF` (`VCF.pm:356-422`) and the records it yields;
    /// the walk ends at a record whose `get_start` lies past the feature's end
    /// plus one (`File.pm:333-335`).
    fn match_record(&self, record: &TabixRecord, group: &RecordGroup<'_>, scores: bool) -> Verdict {
        let spec = self.base.spec();
        // An `exact` source reads a row in whole numbers: its key is the
        // integer that leads the raw `POS` Ensembl VEP hands the allele match
        // (`VCF.pm:416`), and the walk to the row and its name take that key.
        let candidate = if spec.match_type == MatchType::Exact {
            let Some(key) = self.base.keyed_start_of(record) else {
                return Verdict::NoMatch;
            };
            Candidate::of(record, &self.keys, Coordinate::Whole(key), whole_end)
        } else {
            Candidate::of(record, &self.keys, self.base.start_of(record), perl_end)
        };
        if candidate.start.is_past(group.end as i64 + 1) {
            return Verdict::PastEnd;
        }
        let line = &candidate.line;

        // `:364-377`: under `same_type`, or an `exact` match against a
        // breakpoint, the record's structural class is the variant's. The
        // class term of a sequence variant without a stored class is never a
        // breakpoint, so it is computed only when a rule can read it.
        let structural = group.lead.is_structural || group.lead.variant_class.is_structural();
        let class_read = spec.same_type
            || (spec.match_type == MatchType::Exact
                && (structural || group.lead.record_so_term.is_some()));
        if class_read {
            let variant_class = group.lead.class_so_term();
            let exact_breakpoint =
                spec.match_type == MatchType::Exact && variant_class.contains("breakpoint");
            if (spec.same_type || exact_breakpoint)
                && record_so_term(&candidate.alts) != Some(variant_class.as_str())
            {
                return Verdict::NoMatch;
            }
        }

        // `:379-387`: a structural variant or a structural record takes the
        // base predicates on the record's feature coordinates, an insertion's
        // end being its reference's (`File.pm:526-528`).
        if structural || line.routes_to_structural(&candidate.info) {
            let ref_end = if line.alt.contains("INS") {
                candidate.start.plus(line.reference.len() as i64).plus(-1)
            } else {
                candidate.end
            };
            return match base_overlap(
                spec,
                candidate.start,
                ref_end,
                group.start as i64,
                group.end as i64,
            ) {
                Some(percentage) => Verdict::Match(
                    self.create_records(&candidate, None, scores),
                    Some(percentage),
                ),
                None => Verdict::NoMatch,
            };
        }

        match spec.match_type {
            MatchType::Overlap => {
                // `:389-404`: the record spans `POS..POS + length(REF) - 1`
                // and the percentage is the share of the variant it covers.
                let parser_start = line.pos;
                let parser_end = line.reference_end();
                let (vs, ve) = (
                    group.start.min(group.end) as i64,
                    group.start.max(group.end) as i64,
                );
                // The compiled `overlap` reads each end's whole number (`:403`).
                if parser_end.whole() < vs || parser_start.whole() > ve {
                    return Verdict::NoMatch;
                }
                let shared = match (parser_start, parser_end) {
                    (Coordinate::Whole(start), Coordinate::Whole(end)) => {
                        (1 + ve.min(end) - vs.max(start)) as f64
                    }
                    (start, end) => {
                        1.0 + (ve as f64).min(end.value()) - (vs as f64).max(start.value())
                    }
                };
                let percentage = format!("{:.3}", 100.0 * shared / (ve - vs + 1) as f64);
                Verdict::Match(
                    self.create_records(&candidate, None, scores),
                    Some(percentage),
                )
            }
            MatchType::Exact => {
                let matches = exact_matches(line, &candidate.alts, group);
                if matches.is_empty() {
                    return Verdict::NoMatch;
                }
                Verdict::Match(
                    self.create_records(&candidate, Some(&matches), scores),
                    None,
                )
            }
            MatchType::Within | MatchType::Surrounding => Verdict::NoMatch,
        }
    }
}

/// One line of the file by its VCF columns, an absent column reading as
/// empty.
struct VcfLine<'a> {
    chr: &'a str,
    /// `POS` as Perl reads the column or, for an `exact` source, the integer
    /// that leads it, the source's key.
    pos: Coordinate,
    id: &'a str,
    reference: &'a str,
    alt: &'a str,
    qual: &'a str,
    filter: &'a str,
    info: &'a str,
}

impl<'a> VcfLine<'a> {
    fn of(record: &'a TabixRecord) -> Self {
        let column = |index: usize| record.columns.get(index).map_or("", String::as_str);
        VcfLine {
            chr: &record.chr,
            pos: whole_number(record.start),
            id: column(2),
            reference: column(3),
            alt: column(4),
            qual: column(5),
            filter: column(6),
            info: column(7),
        }
    }

    /// `get_alternatives` (`BaseVCF4.pm:486-490`): the ALT column split at
    /// `,`; none for a column Perl reads as false.
    fn alternatives(&self) -> Vec<&'a str> {
        if truthy(self.alt) {
            perl_split(self.alt, ',')
        } else {
            Vec::new()
        }
    }

    /// The INFO values of [`InfoValues`]: a key is an item of the column, the
    /// text between its `;` (or the column's start or end) and its `=` (or the
    /// item's end), so each key is searched for in its delimited forms, and
    /// the last item with the key wins. `my ($key, $value) = split('=', $info)`
    /// (`BaseVCF4.pm:582`) takes the text between the item's first and second
    /// `=`, so `END=` is the defined empty value and a key alone has none.
    fn info_values(&self, keys: &InfoKeys) -> InfoValues<'a> {
        // `split(';', ($self->get_raw_info || ''))`: a column Perl reads as
        // false has no item.
        let text = if truthy(self.info) { self.info } else { "" };
        let info = text.as_bytes();
        // The value of the item whose key ends at `after`.
        let value_at = |after: usize| -> Option<&'a str> {
            if info.get(after) != Some(&b'=') {
                return None;
            }
            let item_end = info[after..]
                .find_byte(b';')
                .map_or(info.len(), |end| after + end);
            Some(text[after + 1..item_end].split('=').next().unwrap_or(""))
        };
        // Whether `key`, ending at `after`, is an item's: the item ends there
        // or its value follows. Perl's `split` drops the trailing empty items,
        // so the empty key names one only when a non-empty item follows it.
        let closes = |key: &[u8], after: usize| {
            (after == info.len() || matches!(info[after], b';' | b'='))
                && (!key.is_empty() || info[after..].iter().any(|&b| b != b';'))
        };
        let value_of = |needle: &[u8]| -> Option<Option<&'a str>> {
            let key = needle.strip_prefix(b";")?;
            let mut last: Option<usize> = None;
            if info.starts_with(key) && closes(key, key.len()) {
                last = Some(key.len());
            }
            // Inside the column every item opens after a `;`.
            for at in info.find_iter(needle) {
                let after = at + needle.len();
                if closes(key, after) {
                    last = Some(after);
                }
            }
            last.map(value_at)
        };
        InfoValues {
            svlen: value_of(&keys.svlen),
            end: value_of(&keys.end),
            svtype: value_of(&keys.svtype),
            fields: keys.fields.iter().map(|needle| value_of(needle)).collect(),
        }
    }

    /// Whether ensembl-io places the record as a structural variant
    /// (`BaseVCF4.pm:262, 310`): `SVTYPE` anywhere in the INFO text, or a
    /// bracket in the ALT column.
    fn structural_shape(&self) -> bool {
        (truthy(self.info) && self.info.contains("SVTYPE"))
            || (truthy(self.alt) && self.alt.contains(['<', '[', ']', '>']))
    }

    /// `get_start` (`BaseVCF4.pm:255-276`): `POS`, one more for a record of
    /// structural shape ([`Self::structural_shape`]) or when one of `alts`,
    /// the record's alternates, differs in length from the reference.
    fn feature_start(&self, alts: &[&str], structural: bool) -> Coordinate {
        let indel = alts.iter().any(|alt| alt.len() != self.reference.len());
        if structural || indel {
            self.pos.plus(1)
        } else {
            self.pos
        }
    }

    /// `get_end` (`BaseVCF4.pm:297-317`): the start plus the first `SVLEN`'s
    /// magnitude less one, else `END` as `end_of` reads it, else the start of
    /// a structural record, else `POS` plus the reference's length less one.
    fn feature_end(
        &self,
        start: Coordinate,
        info: &InfoValues<'_>,
        structural: bool,
        end_of: fn(&str) -> Coordinate,
    ) -> Coordinate {
        if let Some(Some(svlen)) = info.svlen {
            // `(split(',', $svlen))[0]`: the text before the first `,`.
            let first = svlen.split(',').next().unwrap_or("");
            // `as i64` saturates a value past the integer range.
            let length = perl_number(first).abs() as i64;
            return match start {
                Coordinate::Whole(start) => {
                    Coordinate::Whole(start.saturating_add(length).saturating_sub(1))
                }
                float => float.plus(length).plus(-1),
            };
        }
        if let Some(Some(end)) = info.end {
            return end_of(end);
        }
        if structural {
            return start;
        }
        self.reference_end()
    }

    /// `POS` plus the reference's length less one (`BaseVCF4.pm:314`,
    /// `VCF.pm:391`), which holds at the lowest integer for an empty reference
    /// at a `POS` of that integer.
    fn reference_end(&self) -> Coordinate {
        self.pos.plus(self.reference.len() as i64).plus(-1)
    }

    /// `VCF.pm:384-385`: a record whose ALT is not bases alone and that
    /// carries a true `SVTYPE` or a symbolic or breakend allele form.
    fn routes_to_structural(&self, info: &InfoValues<'_>) -> bool {
        let bases = !self.alt.is_empty()
            && self
                .alt
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'));
        if bases {
            return false;
        }
        let svtype = info.svtype.flatten().is_some_and(truthy);
        svtype || breakend_or_symbolic(self.alt)
    }

    /// `_get_score` under `$get_scores and $self->_get_score`
    /// (`VCF.pm:237, 249, 314-318`): QUAL unless `.` or Perl's false.
    fn score(&self) -> Option<CustomValue> {
        (self.qual != "." && truthy(self.qual)).then(|| CustomValue::Text(self.qual.to_string()))
    }
}

/// An `END` as Perl reads it: a plain whole number as it stands, any other
/// cell the number Perl reads from it, the whole number of its float
/// ([`whole_end`]) for an infinity or a text that is no number.
fn perl_end(end: &str) -> Coordinate {
    match end.parse() {
        Ok(whole) => Coordinate::Whole(whole),
        Err(_) => Coordinate::read(end).unwrap_or_else(|| whole_end(end)),
    }
}

/// An `END` in whole numbers, as an `exact` source reads it: the whole number
/// of the float Perl reads from the cell, held at the ends of the 64-bit
/// range.
fn whole_end(end: &str) -> Coordinate {
    Coordinate::Whole(perl_number(end) as i64)
}

/// `/[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$/` (`VCF.pm:385`): a bracketed form
/// with no `*` inside it, or a single breakend.
fn breakend_or_symbolic(alt: &str) -> bool {
    let bytes = alt.as_bytes();
    let bracketed = (0..bytes.len()).any(|i| {
        matches!(bytes[i], b'<' | b'[' | b']')
            && bytes[i + 1..]
                .iter()
                .take_while(|&&b| b != b'*')
                .enumerate()
                .any(|(k, &b)| k >= 1 && matches!(b, b'>' | b']' | b'['))
    });
    let word = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    let opens_breakend = bytes.len() >= 2 && bytes[0] == b'.' && word(bytes[1]);
    let closes_breakend =
        bytes.len() >= 2 && bytes[bytes.len() - 1] == b'.' && word(bytes[bytes.len() - 2]);
    bracketed || opens_breakend || closes_breakend
}

/// `get_SO_term` (`Parser.pm:663-728`) on the record's alternates joined by
/// `,`: the SO term of its structural type, `None` for a sequence record or
/// a type the term table lacks.
fn record_so_term(alts: &[&str]) -> Option<&'static str> {
    let joined = alts.join(",");
    sv_so_term_and_class(&sv_type_abbreviation(&joined, None)).map(|(term, _)| term)
}

/// `get_matched_variant_alleles` (`Utils/Sequence.pm:1144-1245`) between the
/// variant's alternates and the record's: each pair of reference and
/// alternate is trimmed to its minimal form from the left then the right and,
/// when either allele is longer than one base, from the right then the left
/// too; an alternate of the record matches an alternate of the variant when a
/// trimmed form and its start agree. The variant's alleles are read on the
/// forward strand (`:1182, 1196`). Each matched record alternate is returned
/// with its index in the ALT column and the variant it matched, in ALT order.
fn exact_matches<'a>(
    line: &VcfLine<'_>,
    alts: &[&str],
    group: &RecordGroup<'a>,
) -> Vec<(usize, &'a InputVariant)> {
    let pos = line.pos.whole();
    if pos == 0 {
        return Vec::new();
    }
    let mut by_minimal_form: HashMap<(Vec<u8>, Vec<u8>, i64), &'a InputVariant> = HashMap::new();
    for variant in &group.alleles {
        if variant.start == 0 {
            continue;
        }
        let (reference, alt) = variant.forward_strand_alleles();
        let (reference, alt) = (dashed(&reference), dashed(&alt));
        for end_first in trim_directions(&reference, &alt) {
            let form = trim_sequences(&reference, &alt, variant.start as i64, *end_first);
            by_minimal_form.insert(form, variant);
        }
    }
    let mut matches = Vec::new();
    for (index, alt) in alts.iter().enumerate() {
        let (reference, alt) = (line.reference.as_bytes(), alt.as_bytes());
        for end_first in trim_directions(reference, alt) {
            let form = trim_sequences(reference, alt, pos, *end_first);
            if let Some(variant) = by_minimal_form.get(&form) {
                matches.push((index, *variant));
                break;
            }
        }
    }
    matches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::custom::grammar::parse_custom;
    use crate::custom::source::CustomSources;
    use crate::test_fixtures::{write_bgzf_bytes, write_tabix_fixture, IndexSpec};
    use vep_core::coordinate::Strand;
    use vep_core::custom_annotation::{CustomAnnotations, CustomStats};
    use vep_core::variant::VariantClass;

    const HEADER: [&str; 2] = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ];

    fn fixture(dir: &Path, name: &str, header: &[&str], records: &[&str]) -> String {
        let lines: Vec<&str> = header.iter().chain(records).copied().collect();
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    fn sources(argument: &str, options: CustomRunOptions) -> CustomSources {
        sources_of(&[argument], options)
    }

    /// The sources of several `--custom` arguments, in the order given.
    fn sources_of(arguments: &[&str], options: CustomRunOptions) -> CustomSources {
        CustomSources::open(
            arguments
                .iter()
                .map(|argument| parse_custom(argument).unwrap())
                .collect(),
            Arc::default(),
            options,
        )
        .unwrap()
    }

    /// A chopped allele of one input record, as the parser lays the record's
    /// alleles out: the same span and record ordinal for every allele.
    fn allele(start: u64, end: u64, reference: &str, alt: &str) -> InputVariant {
        let mut variant = InputVariant::new(
            "21".into(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.input_record = 1;
        variant
    }

    fn annotate(sources: &CustomSources, batch: &mut [InputVariant]) {
        sources.run_batch(batch).unwrap();
    }

    /// Each variant as the only allele of its own input record.
    fn records(mut batch: Vec<InputVariant>) -> Vec<InputVariant> {
        for (i, variant) in batch.iter_mut().enumerate() {
            variant.input_record = i as u64 + 1;
        }
        batch
    }

    fn only(variant: &InputVariant) -> &CustomAnnotations {
        assert_eq!(variant.custom_annotations.len(), 1, "{variant:?}");
        &variant.custom_annotations[0]
    }

    fn text(field: &CustomFieldValue) -> &str {
        match field {
            CustomFieldValue::Text(text) => text,
            other => panic!("{other:?} is not a text"),
        }
    }

    /// A listed record as its name and its text fields other than `FILTER`.
    type Named = (String, Vec<(String, String)>);

    /// The records every entry of `variant` lists, in entry order.
    fn listed(variant: &InputVariant) -> Vec<Named> {
        variant
            .custom_annotations
            .iter()
            .flat_map(|entry| entry.records.iter())
            .map(|r| {
                (
                    r.name.clone(),
                    r.fields
                        .iter()
                        .filter(|(k, _)| *k != "FILTER")
                        .map(|(k, v)| (k.clone(), text(v).to_string()))
                        .collect(),
                )
            })
            .collect()
    }

    fn named(name: &str, fields: &[(&str, &str)]) -> Named {
        (
            name.to_string(),
            fields
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect(),
        )
    }

    /// `BaseVCF4.pm:576-589` through `info_values`: a key's value is the text
    /// between its first and second `=`, so a key ending in `=` has the empty
    /// value, a key alone has none, a later item replaces an earlier, a key
    /// inside another key's text is no item, and a key the column lacks is
    /// absent; `SVLEN`, `END` and `SVTYPE` are read the same way, as fields
    /// too when the source names them. The empty key (`fields=%A`) names the
    /// empty items Perl's `split` keeps (`A;;B`, `;A`) and not the trailing
    /// ones it drops, and a column Perl reads as false has no item at all.
    #[test]
    fn info_values_follow_perls_split() {
        let record = |reference: &str, alt: &str, info: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: ["", "100", ".", reference, alt, "50", "PASS", info]
                .iter()
                .map(|c| c.to_string())
                .collect(),
            header: None,
        };
        let plain = record("A", "G", "DB;AF=0.1;X=;Y==5;AF=0.2;Z=a=b");
        let line = VcfLine::of(&plain);
        let fields: Vec<String> = ["DB", "AF", "X", "Y", "Z", "Q"]
            .iter()
            .map(|f| f.to_string())
            .collect();
        let info = line.info_values(&InfoKeys::new(&fields));
        assert_eq!(
            info.fields,
            [
                Some(None),
                Some(Some("0.2")),
                Some(Some("")),
                Some(Some("")),
                Some(Some("a")),
                None
            ]
        );
        assert_eq!((info.svlen, info.end, info.svtype), (None, None, None));
        assert!(!line.structural_shape());
        assert_eq!(line.feature_start(&line.alternatives(), false).whole(), 100);
        assert_eq!(
            line.feature_end(100.into(), &info, false, perl_end).whole(),
            100
        );
        let structural = record("N", "<DEL>", "SVTYPE=DEL;END=150;SVLEN=-50");
        let line = VcfLine::of(&structural);
        let keys: Vec<String> = ["END", "SVTYPE", "XEND", "END;SVTYPE"]
            .iter()
            .map(|f| f.to_string())
            .collect();
        let info = line.info_values(&InfoKeys::new(&keys));
        assert_eq!(info.svtype, Some(Some("DEL")));
        assert_eq!(
            info.fields,
            [Some(Some("150")), Some(Some("DEL")), None, None]
        );
        assert!(line.structural_shape());
        assert_eq!(
            line.feature_end(101.into(), &info, true, perl_end).whole(),
            150
        );
        assert!(line.routes_to_structural(&info));

        // `END` against a column, and the fields `` and `0` against it.
        let read = |info: &str| -> (Option<Option<String>>, Vec<Option<Option<String>>>) {
            let record = record("A", "G", info);
            let line = VcfLine::of(&record);
            let values = line.info_values(&InfoKeys::new(&[String::new(), "0".to_string()]));
            let owned = |value: Option<Option<&str>>| value.map(|v| v.map(str::to_string));
            (
                owned(values.end),
                values.fields.iter().map(|&v| owned(v)).collect(),
            )
        };
        assert_eq!(read("XEND=1;ENDX=2;END=3").0, Some(Some("3".to_string())));
        assert_eq!(read("XEND=1").0, None);
        assert_eq!(read("A;;B").1, [Some(None), None]);
        assert_eq!(read(";A").1, [Some(None), None]);
        assert_eq!(read("=5").1, [Some(Some("5".to_string())), None]);
        assert_eq!(read("0=1").1, [None, Some(Some("1".to_string()))]);
        for column in ["", "A;", "A;;", ";", "0"] {
            assert_eq!(read(column).1, [None, None], "{column:?}");
        }
    }

    /// `VCF.pm:364-377` under `same_type=1`: a record matches only when
    /// `get_SO_term` of its alternates is the variant's class term, so a
    /// `<DEL>` record matches a `<DEL>` variant and not a `<DUP>` one, and a
    /// sequence record, whose term is undefined (`Parser.pm:722`), matches no
    /// variant at all; an `exact` source matches a breakpoint variant only
    /// against a breakend record (`:369`).
    #[test]
    fn same_type_admits_structural_records_of_the_variants_class_alone() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "types.vcf.gz",
            &HEADER,
            &[
                "21\t100\tsnv\tA\tG\t.\tPASS\tDP=1",
                "21\t300\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=400",
                "21\t500\tbnd\tN\tN[22:1000[\t.\tPASS\tSVTYPE=BND",
            ],
        );
        let typed = sources(
            &format!("file={file},short_name=S,format=vcf,type=overlap,same_type=1"),
            CustomRunOptions::default(),
        );
        // A structural record's variant as the parser builds it: its class and
        // the record's class term set.
        let structural = |start: u64, end: u64, alt: &str, class: VariantClass, record: u64| {
            let mut variant = allele(start, end, "N", alt);
            variant.is_structural = true;
            variant.variant_class = class;
            variant.record_so_term = Some(class.so_term().to_string());
            variant.input_record = record;
            variant
        };
        let mut batch = vec![
            allele(100, 100, "A", "G"),
            structural(301, 400, "<DEL>", VariantClass::StructuralDeletion, 2),
            structural(301, 400, "<DUP>", VariantClass::Duplication, 3),
        ];
        annotate(&typed, &mut batch);
        assert!(batch[0].custom_annotations.is_empty());
        assert_eq!(only(&batch[1]).records[0].name, "del");
        assert!(batch[2].custom_annotations.is_empty());

        let exact = sources(
            &format!("file={file},short_name=S,format=vcf,type=exact"),
            CustomRunOptions::default(),
        );
        let mut batch = vec![
            structural(501, 501, "N[22:1000[", VariantClass::Translocation, 1),
            structural(301, 400, "<DEL>", VariantClass::StructuralDeletion, 2),
        ];
        annotate(&exact, &mut batch);
        assert_eq!(only(&batch[0]).records[0].name, "bnd");
        assert_eq!(only(&batch[1]).records[0].name, "del");
    }

    /// `Utils/Sequence.pm:1199-1240`: the record `A G,T` at 100 matches the
    /// variant `A/T` through its second alternate, index 1, whose per-allele
    /// value the record takes (`VCF.pm:241`), and the variant's own allele
    /// is the record's `allele` (`:235`); the record `CA CG` at 99 matches
    /// `A/G` at 100 once its shared `C` is trimmed and its start moved up;
    /// `TA T` at 99 matches the deletion `A/-` at 100 the same way; a variant
    /// given on the reverse strand is compared through its forward-strand
    /// alleles (`:1182, 1196`); an alternate matching no allele yields no
    /// record.
    #[test]
    fn exact_matches_each_alt_by_trimmed_alleles_and_returns_its_index() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "exact.vcf.gz",
            &HEADER,
            &[
                "21\t99\tanchored\tCA\tCG\t.\tPASS\tAF=0.3",
                "21\t99\tdel\tTA\tT\t.\tPASS\tAF=0.4",
                "21\t100\tmulti\tA\tG,T\t.\tPASS\tAF=0.1,0.2",
                "21\t100\tother\tA\tC\t.\tPASS\tAF=0.5",
            ],
        );
        let sources = sources(
            &format!("file={file},short_name=S,format=vcf,type=exact,fields=AF"),
            CustomRunOptions::default(),
        );

        let mut batch = vec![allele(100, 100, "A", "T")];
        annotate(&sources, &mut batch);
        let records = &only(&batch[0]).records;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "multi");
        assert_eq!(records[0].allele.as_deref(), Some(&b"T"[..]));
        assert_eq!(text(&records[0].fields["AF"]), "0.2");

        let mut batch = vec![allele(100, 100, "A", "G")];
        annotate(&sources, &mut batch);
        let names: Vec<(&str, &str)> = only(&batch[0])
            .records
            .iter()
            .map(|r| (r.name.as_str(), text(&r.fields["AF"])))
            .collect();
        assert_eq!(names, [("anchored", "0.3"), ("multi", "0.1")]);

        let mut batch = vec![allele(100, 100, "A", "-")];
        annotate(&sources, &mut batch);
        assert_eq!(only(&batch[0]).records[0].name, "del");

        let mut reverse = allele(100, 100, "T", "A");
        reverse.strand = Strand::Reverse;
        let mut batch = vec![reverse];
        annotate(&sources, &mut batch);
        assert_eq!(only(&batch[0]).records[0].name, "multi");
        assert_eq!(
            only(&batch[0]).records[0].allele.as_deref(),
            Some(&b"A"[..])
        );

        let mut batch = vec![allele(100, 100, "A", "N")];
        annotate(&sources, &mut batch);
        assert!(batch[0].custom_annotations.is_empty());
    }

    /// `VCF.pm:389-422` has branches for `overlap` and `exact` alone, so a
    /// sequence record never matches a `within` or `surrounding` source,
    /// however it lies against the variant; the same record matches an
    /// `overlap` source over `POS..POS + length(REF) - 1` with the share of
    /// the variant it covers (`:401`), and a structural record takes the base
    /// predicates (`:379-387`).
    #[test]
    fn small_variant_within_and_surrounding_never_match() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "small.vcf.gz",
            &HEADER,
            &[
                "21\t100\tsnv\tA\tG\t.\tPASS\tDP=5",
                "21\t150\tdel\tCATG\tC\t.\tPASS\tDP=6",
                "21\t300\tsv\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=400",
            ],
        );
        for kind in ["within", "surrounding"] {
            let sources = sources(
                &format!("file={file},short_name=S,format=vcf,type={kind},fields=DP"),
                CustomRunOptions::default(),
            );
            let mut batch = records(vec![
                allele(100, 100, "A", "C"),
                allele(140, 160, "ACATGAAAAAAAAAAAAAAAA", "A"),
                allele(151, 152, "AT", "-"),
            ]);
            annotate(&sources, &mut batch);
            assert!(
                batch.iter().all(|v| v.custom_annotations.is_empty()),
                "{kind}"
            );
            let mut structural = allele(301, 400, "N", "<DEL>");
            structural.is_structural = true;
            structural.input_record = 2;
            let mut batch = vec![structural];
            annotate(&sources, &mut batch);
            assert_eq!(only(&batch[0]).records[0].name, "sv", "{kind}");
        }
        let overlap = sources(
            &format!("file={file},short_name=S,format=vcf,type=overlap,fields=DP%PC"),
            CustomRunOptions::default(),
        );
        let mut batch = records(vec![
            allele(100, 100, "A", "C"),
            allele(152, 153, "TG", "-"),
            allele(154, 154, "A", "T"),
        ]);
        annotate(&overlap, &mut batch);
        assert_eq!(only(&batch[0]).records[0].name, "snv");
        assert_eq!(text(&only(&batch[0]).records[0].fields["PC"]), "100.000");
        assert_eq!(only(&batch[1]).records[0].name, "del");
        assert_eq!(text(&only(&batch[1]).records[0].fields["PC"]), "100.000");
        assert!(batch[2].custom_annotations.is_empty());
    }

    /// `File.pm:526-528`: a structural VCF record whose ALT names `INS` ends
    /// at its start plus the reference's length less one, whatever `END`
    /// says, so a variant inside the `END` span but past the reference does
    /// not overlap it; the record's name and the walk still read ensembl-io's
    /// coordinates (`BaseVCF4.pm:255-317`): a deletion's name opens one past
    /// `POS`, and `SVLEN` sets the end before `END` does (`:302-308`).
    #[test]
    fn ins_records_end_at_pos_plus_ref_length() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "sv.vcf.gz",
            &HEADER,
            &[
                "21\t100\t.\tA\t<INS>\t.\tPASS\tSVTYPE=INS;END=150",
                "21\t200\t.\tAT\tA\t.\tPASS\tDP=1",
                "21\t300\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=360",
            ],
        );
        let sources = sources(
            &format!("file={file},format=vcf,type=overlap,coords=1"),
            CustomRunOptions::default(),
        );
        let mut batch = records(vec![
            allele(101, 101, "A", "G"),
            allele(120, 120, "A", "G"),
            allele(201, 201, "T", "-"),
            allele(340, 340, "A", "G"),
            allele(355, 355, "A", "G"),
        ]);
        annotate(&sources, &mut batch);
        assert_eq!(only(&batch[0]).records[0].name, "21:101-150");
        assert!(batch[1].custom_annotations.is_empty());
        assert_eq!(only(&batch[2]).records[0].name, "21:201-201");
        assert_eq!(only(&batch[3]).records[0].name, "21:301-350");
        assert!(batch[4].custom_annotations.is_empty());
    }

    /// `get_end` on `SVLEN` (`BaseVCF4.pm:302-305`) is the start plus the
    /// length's magnitude less one. From a whole start the sum is held at the
    /// upper end of the 64-bit range before the one comes off: a length of
    /// 9223372036854775807 from a start of 101 ends at 9223372036854775806.
    /// A start with a fraction adds in the float: a deletion at `POS` 100.5
    /// with `SVLEN=-50` spans 101.5 to 150.5, which Ensembl VEP 116.2 lists
    /// on the pinned image from 101 to 150 for the record at base 101, 1
    /// percent of it covered, and not for the record at base 100, whose end
    /// plus one its start lies past (`File.pm:333-335`).
    #[test]
    fn an_svlen_end_saturates_from_a_whole_start_and_keeps_the_fraction_of_a_float_start() {
        let dir = tempfile::tempdir().unwrap();
        let listed_from = |name: &str, row: &str| -> Vec<Vec<Named>> {
            let file = fixture(dir.path(), name, &HEADER, &[row]);
            let sources = sources(
                &format!("file={file},format=vcf,type=overlap,coords=1,fields=PC"),
                CustomRunOptions::default(),
            );
            let mut batch = records(vec![allele(100, 100, "A", "G"), allele(101, 101, "A", "G")]);
            annotate(&sources, &mut batch);
            batch.iter().map(listed).collect()
        };
        assert_eq!(
            listed_from(
                "whole.vcf.gz",
                "21\t100\t.\tN\t<DEL>\t.\t.\tSVLEN=9223372036854775807"
            ),
            [
                vec![],
                vec![named("21:101-9223372036854775806", &[("PC", "0.000")])]
            ]
        );
        assert_eq!(
            listed_from("fraction.vcf.gz", "21\t100.5\t.\tN\t<DEL>\t.\t.\tSVLEN=-50"),
            [vec![], vec![named("21:101-150", &[("PC", "1.000")])]]
        );
    }

    /// An `END` Perl reads as an infinity or as no number is no coordinate
    /// ([`Coordinate::read`]), so the row ends at the whole number of that
    /// float: the record at the start of a deletion at `POS` 100 lists it as
    /// `21:101-9223372036854775807` for `END=inf`, the upper end of the 64-bit
    /// range, and lists nothing for `END=nan`, which ends it at 0. Ensembl VEP
    /// 116.2 lists neither row, its reader returning each for no query
    /// (`tabix_row_with_an_end_below_its_start`).
    #[test]
    fn an_end_perl_reads_as_an_infinity_or_as_no_number_is_the_whole_number_of_its_float() {
        let dir = tempfile::tempdir().unwrap();
        let listed_from = |name: &str, end: &str| -> Vec<Vec<Named>> {
            let row = format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END={end}");
            let file = fixture(dir.path(), name, &HEADER, &[row.as_str()]);
            let sources = sources(
                &format!("file={file},format=vcf,type=overlap,coords=1"),
                CustomRunOptions::default(),
            );
            let mut batch = records(vec![allele(100, 100, "A", "G"), allele(101, 101, "A", "G")]);
            annotate(&sources, &mut batch);
            batch.iter().map(listed).collect()
        };
        assert_eq!(
            listed_from("inf.vcf.gz", "inf"),
            [vec![], vec![named("21:101-9223372036854775807", &[])]]
        );
        assert_eq!(listed_from("nan.vcf.gz", "nan"), [vec![], vec![]]);
    }

    /// A row whose INFO `END` lies below its position, or is not a number
    /// (`END=x`, `END=.`, `END=`), spans its reference allele, as htslib 1.11
    /// and later span it, so a variant inside that allele is matched by
    /// `overlap` and the row's own deletion by `exact`; the two well-formed
    /// rows are the controls. Ensembl VEP 116.2's reader (htslib 1.9,
    /// `tbx.c:133`, `hts.c:2645`) gives such a row a span that ends before it
    /// starts and returns it for no query over its bases
    /// (`tabix_row_with_an_end_below_its_start`).
    #[test]
    fn a_row_whose_end_lies_below_its_position_spans_its_reference_allele() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "span.vcf.gz",
            &HEADER,
            &[
                "21\t5029100\tctrl_no_end\tTCTG\tT\t.\t.\tDP=3",
                "21\t5029120\tctrl_end\tATGA\tA\t.\t.\tEND=5029123",
                "21\t5029140\tend_below\tCAGC\tC\t.\t.\tEND=5029000",
                "21\t5029160\tend_x\tCCAG\tC\t.\t.\tEND=x",
                "21\t5029200\tend_dot\tAGAA\tA\t.\t.\tEND=.",
                "21\t5029220\tend_empty\tTCGT\tT\t.\t.\tEND=",
            ],
        );
        let names = [
            "ctrl_no_end",
            "ctrl_end",
            "end_below",
            "end_x",
            "end_dot",
            "end_empty",
        ];
        let first_record_names = |batch: &[InputVariant]| -> Vec<String> {
            batch
                .iter()
                .map(|variant| only(variant).records[0].name.clone())
                .collect()
        };
        let overlap = sources(
            &format!("file={file},short_name=SPAN,format=vcf,type=overlap,coords=0"),
            CustomRunOptions::default(),
        );
        let mut substitutions = records(vec![
            allele(5029102, 5029102, "T", "G"),
            allele(5029122, 5029122, "G", "T"),
            allele(5029142, 5029142, "G", "T"),
            allele(5029162, 5029162, "A", "C"),
            allele(5029202, 5029202, "A", "C"),
            allele(5029222, 5029222, "G", "T"),
        ]);
        annotate(&overlap, &mut substitutions);
        assert_eq!(first_record_names(&substitutions), names);
        let exact = sources(
            &format!("file={file},short_name=SPAN,format=vcf,type=exact,coords=0"),
            CustomRunOptions::default(),
        );
        let mut deletions = records(vec![
            allele(5029101, 5029103, "CTG", "-"),
            allele(5029121, 5029123, "TGA", "-"),
            allele(5029141, 5029143, "AGC", "-"),
            allele(5029161, 5029163, "CAG", "-"),
            allele(5029201, 5029203, "GAA", "-"),
            allele(5029221, 5029223, "CGT", "-"),
        ]);
        annotate(&exact, &mut deletions);
        assert_eq!(first_record_names(&deletions), names);
    }

    /// A row whose `POS` is a minus sign and digits (`-5`) is matched as
    /// Ensembl VEP 116.2 matches it. `tabix` indexes it from the first base
    /// over the length of its reference allele (htslib 1.9 `tbx.c:99`, `:123`)
    /// and Ensembl VEP reads the `POS` written (`BaseVCF4.pm:244-247`): under
    /// `overlap` a row at `-5` whose reference allele is ten bases spans -5 to
    /// 4 (`VCF.pm:389-404`), so the records at bases 1 and 4 list it, under
    /// its ID or, with `coords`, as `21:-4-4` (`BaseVCF4.pm:255-317`), and the
    /// record at base 5 does not; a row of one base there reaches no record,
    /// and under `exact` no record starts where either row does.
    #[test]
    fn a_row_whose_pos_a_minus_sign_leads_is_matched_on_the_pos_written() {
        let dir = tempfile::tempdir().unwrap();
        let ten_bases = "A".repeat(10);
        let below_zero = format!("21\t-5\tbelow\t{ten_bases}\tA\t.\t.\tAF=0.2");
        let file = fixture(
            dir.path(),
            "below.vcf.gz",
            &HEADER,
            &[
                "21\t-5\tone_base\tA\tC\t.\t.\tAF=0.1",
                &below_zero,
                "21\t100\tctrl\tA\tG\t.\t.\tAF=0.3",
            ],
        );
        let listed_under = |options: &str| -> Vec<Vec<Named>> {
            let sources = sources(
                &format!("file={file},format=vcf,{options}"),
                CustomRunOptions::default(),
            );
            let mut batch = records(vec![
                allele(1, 1, "A", "G"),
                allele(4, 4, "A", "G"),
                allele(5, 5, "A", "G"),
                allele(100, 100, "A", "G"),
            ]);
            annotate(&sources, &mut batch);
            batch.iter().map(listed).collect()
        };
        let frequency = |name: &str, frequency: &str| vec![named(name, &[("AF", frequency)])];
        assert_eq!(
            listed_under("type=overlap,coords=0,fields=AF"),
            [
                frequency("below", "0.2"),
                frequency("below", "0.2"),
                vec![],
                frequency("ctrl", "0.3"),
            ]
        );
        let whole = |name: &str, frequency: &str| {
            vec![named(name, &[("AF", frequency), ("PC", "100.000")])]
        };
        assert_eq!(
            listed_under("type=overlap,coords=1,fields=AF%PC"),
            [
                whole("21:-4-4", "0.2"),
                whole("21:-4-4", "0.2"),
                vec![],
                whole("21:100-100", "0.3"),
            ]
        );
        assert_eq!(
            listed_under("type=exact,coords=0,fields=AF"),
            [vec![], vec![], vec![], frequency("ctrl", "0.3")]
        );
    }

    /// Under `exact` a row is keyed by the leading integer of its `POS`
    /// ([`TabixSource::keyed_start_of`]): a row at `0.5e1`, which the reader
    /// places at base 0 and Perl reads as 5, matches no record, so the record
    /// at base 5 that its trimmed alleles spell lists the row at `5` alone.
    #[test]
    fn exact_keys_a_row_by_the_leading_integer_of_its_pos() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "leading.vcf.gz",
            &HEADER,
            &[
                "21\t0.5e1\texponent\tACCCC\tGCCCC\t.\t.\tAF=0.4",
                "21\t5\tplain\tA\tG\t.\t.\tAF=0.3",
            ],
        );
        let sources = sources(
            &format!("file={file},format=vcf,type=exact,coords=0,fields=AF"),
            CustomRunOptions::default(),
        );
        let mut batch = vec![allele(5, 5, "A", "G")];
        annotate(&sources, &mut batch);
        assert_eq!(listed(&batch[0]), [named("plain", &[("AF", "0.3")])]);
    }

    /// An `exact` source reads a row in whole numbers from the integer that
    /// leads its `POS`: rows at `1e2` and `2000e-1` are listed for the records
    /// at bases 1 and 2000 and named by those integers, rows at `+100` and
    /// `100.5` for the record at base 100, and the record at base 200, the
    /// number Perl reads from `2000e-1`, lists nothing. A deletion row at
    /// `1e2` whose `END` is `300.5` is the deletion of bases 2 to 300.
    #[test]
    fn an_exact_source_reads_a_row_by_the_integer_that_leads_its_pos() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "keyed.vcf.gz",
            &HEADER,
            &[
                "21\t1\tplain1\tA\tG\t.\t.\tAF=0.1",
                "21\t1e2\texp\tA\tG\t.\t.\tAF=0.2",
                "21\t1e2\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=300.5",
                "21\t100\tplain100\tA\tG\t.\t.\tAF=0.3",
                "21\t+100\tplus\tA\tG\t.\t.\tAF=0.4",
                "21\t100.5\tfrac\tA\tG\t.\t.\tAF=0.5",
                "21\t2000\tplain2000\tA\tG\t.\t.\tAF=0.6",
                "21\t2000e-1\tnegexp\tA\tG\t.\t.\tAF=0.7",
            ],
        );
        let names_under = |coords: &str| -> Vec<Vec<String>> {
            let sources = sources(
                &format!("file={file},format=vcf,type=exact,coords={coords}"),
                CustomRunOptions::default(),
            );
            let mut deletion = allele(2, 300, "N", "<DEL>");
            deletion.is_structural = true;
            deletion.variant_class = VariantClass::StructuralDeletion;
            deletion.record_so_term = Some(deletion.variant_class.so_term().to_string());
            let mut batch = records(vec![
                allele(1, 1, "A", "G"),
                allele(100, 100, "A", "G"),
                allele(200, 200, "A", "G"),
                allele(2000, 2000, "A", "G"),
                deletion,
            ]);
            annotate(&sources, &mut batch);
            batch
                .iter()
                .map(|variant| listed(variant).into_iter().map(|(name, _)| name).collect())
                .collect()
        };
        assert_eq!(
            names_under("0"),
            [
                vec!["plain1", "exp"],
                vec!["plain100", "plus", "frac"],
                vec![],
                vec!["plain2000", "negexp"],
                vec!["del"],
            ]
        );
        assert_eq!(
            names_under("1"),
            [
                vec!["21:1-1"; 2],
                vec!["21:100-100"; 3],
                vec![],
                vec!["21:2000-2000"; 2],
                vec!["21:2-300"],
            ]
        );
    }

    /// A row whose `POS` is the lowest integer and whose reference cell is
    /// empty ends a base below its `POS` (`VCF.pm:391`), where no integer
    /// lies: its end saturates, so under `overlap` the record at base 1,
    /// which the index hands the row to, lists nothing, as Ensembl VEP 116.2
    /// lists nothing for it.
    #[test]
    fn a_row_at_the_lowest_integer_with_no_reference_base_matches_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "lowest.vcf.gz",
            &HEADER,
            &[
                "21\t-9223372036854775808\tno_base\t\tC\t.\t.\tAF=0.1",
                "21\t100\tctrl\tA\tG\t.\t.\tAF=0.3",
            ],
        );
        let sources = sources(
            &format!("file={file},format=vcf,type=overlap,coords=0,fields=AF"),
            CustomRunOptions::default(),
        );
        let mut batch = records(vec![allele(1, 1, "A", "G"), allele(100, 100, "A", "G")]);
        annotate(&sources, &mut batch);
        assert_eq!(
            batch.iter().map(listed).collect::<Vec<_>>(),
            [vec![], vec![named("ctrl", &[("AF", "0.3")])]]
        );
    }

    /// `File.pm:405-415`: a second record named like a listed `COSV` record
    /// lends it the fields it lacks and is not listed again, while it still
    /// counts (`:339`).
    #[test]
    fn cosv_records_merge_by_name() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "cosmic_synthetic.vcf.gz",
            &HEADER,
            &[
                "21\t100\tCOSV1\tA\tG\t.\tPASS\tGENE=X;CDS=c.1A>G",
                "21\t100\tCOSV1\tA\tG\t.\tPASS\tGENE=Y;AA=p.1",
                "21\t100\tCOSV2\tA\tG\t.\tPASS\tGENE=Z",
            ],
        );
        let sources = sources(
            &format!(
                "file={file},short_name=COSMIC,format=vcf,type=exact,fields=GENE%CDS%AA,summary_stats=count"
            ),
            CustomRunOptions::default(),
        );
        let mut batch = vec![allele(100, 100, "A", "G")];
        annotate(&sources, &mut batch);
        let annotation = only(&batch[0]);
        assert_eq!(
            listed(&batch[0]),
            [
                named("COSV1", &[("GENE", "X"), ("CDS", "c.1A>G"), ("AA", "p.1")]),
                named("COSV2", &[("GENE", "Z")]),
            ]
        );
        assert_eq!(
            annotation.stats,
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(3.0))]
                    .into_iter()
                    .collect()
            ))
        );
    }

    /// `File.pm:406` greps the one record list the short name owns, which
    /// every source of that name extends (`:437`): a `COSV` record a later
    /// source reads under the short name lends the namesake an earlier source
    /// listed the fields it lacks and is not listed again (`:405-415`), by
    /// name alone, whichever allele either matched and under either match
    /// type; the earlier source's values stay. The lent fields reach every
    /// allele of the input record.
    #[test]
    fn cosv_namesakes_merge_across_sources_sharing_a_short_name() {
        let dir = tempfile::tempdir().unwrap();
        let a = fixture(
            dir.path(),
            "a.vcf.gz",
            &HEADER,
            &["21\t100\tCOSV1\tA\tC\t.\tPASS\tGENE=X"],
        );
        let b = fixture(
            dir.path(),
            "b.vcf.gz",
            &HEADER,
            &["21\t100\tCOSV1\tA\tC\t.\tq10\tGENE=Y;AA=p.1"],
        );
        let g = fixture(
            dir.path(),
            "g.vcf.gz",
            &HEADER,
            &["21\t100\tCOSV1\tA\tG\t.\tPASS\tGENE=Y;AA=p.1"],
        );
        let cosmic = |file: &str, match_type: &str| {
            format!("file={file},short_name=COSMIC,format=vcf,type={match_type},fields=GENE%AA")
        };
        let run = |arguments: [&str; 2]| -> Vec<InputVariant> {
            let sources = sources_of(&arguments, CustomRunOptions::default());
            let mut batch = vec![allele(100, 100, "A", "C"), allele(100, 100, "A", "G")];
            annotate(&sources, &mut batch);
            assert_eq!(batch[0].custom_annotations, batch[1].custom_annotations);
            batch
        };

        let batch = run([&cosmic(&a, "exact"), &cosmic(&b, "exact")]);
        assert_eq!(batch[0].custom_annotations.len(), 1, "{batch:?}");
        assert_eq!(
            listed(&batch[0]),
            [named("COSV1", &[("GENE", "X"), ("AA", "p.1")])]
        );
        assert_eq!(
            only(&batch[0]).records[0].fields.get("FILTER"),
            Some(&CustomFieldValue::List(vec!["PASS".to_string()]))
        );

        let batch = run([&cosmic(&b, "exact"), &cosmic(&a, "exact")]);
        assert_eq!(
            listed(&batch[0]),
            [named("COSV1", &[("GENE", "Y"), ("AA", "p.1")])]
        );

        let batch = run([&cosmic(&a, "exact"), &cosmic(&g, "exact")]);
        assert_eq!(
            listed(&batch[0]),
            [named("COSV1", &[("GENE", "X"), ("AA", "p.1")])]
        );
        assert_eq!(
            only(&batch[0]).records[0].allele.as_deref(),
            Some(b"C".as_slice())
        );

        let batch = run([&cosmic(&a, "overlap"), &cosmic(&b, "overlap")]);
        assert_eq!(
            listed(&batch[0]),
            [named("COSV1", &[("GENE", "X"), ("AA", "p.1")])]
        );
        assert_eq!(only(&batch[0]).records[0].allele, None);
    }

    /// Two sources under one short name list their records in source order
    /// (`File.pm:437`): `COSV` records of different names are both listed, a
    /// record one source alone holds is listed once with that source's
    /// fields, and namesakes not named `COSV` are listed twice (`:405`).
    #[test]
    fn sources_sharing_a_short_name_list_their_distinct_records_in_order() {
        let dir = tempfile::tempdir().unwrap();
        let file = |name: &str, line: &str| fixture(dir.path(), name, &HEADER, &[line]);
        let a = file("a.vcf.gz", "21\t100\tCOSV1\tA\tC\t.\tPASS\tGENE=X");
        let c = file("c.vcf.gz", "21\t100\tCOSV2\tA\tC\t.\tPASS\tGENE=Z");
        let n = file("n.vcf.gz", "21\t110\tCOSV1\tA\tC\t.\tPASS\tGENE=X");
        let e = file("e.vcf.gz", "21\t100\trsX\tA\tC\t.\tPASS\tGENE=X");
        let f = file("f.vcf.gz", "21\t100\trsX\tA\tC\t.\tPASS\tGENE=Y;AA=p.1");
        let run = |first: &str, second: &str| -> (usize, Vec<Named>) {
            let cosmic = |file: &str| {
                format!("file={file},short_name=COSMIC,format=vcf,type=exact,fields=GENE%AA")
            };
            let sources = sources_of(
                &[&cosmic(first), &cosmic(second)],
                CustomRunOptions::default(),
            );
            let mut batch = vec![allele(100, 100, "A", "C")];
            annotate(&sources, &mut batch);
            (batch[0].custom_annotations.len(), listed(&batch[0]))
        };
        let x = [("GENE", "X")];
        assert_eq!(
            run(&a, &c),
            (
                2,
                vec![named("COSV1", &x), named("COSV2", &[("GENE", "Z")])]
            )
        );
        assert_eq!(
            run(&c, &a),
            (
                2,
                vec![named("COSV2", &[("GENE", "Z")]), named("COSV1", &x)]
            )
        );
        assert_eq!(run(&a, &n), (1, vec![named("COSV1", &x)]));
        assert_eq!(
            run(&e, &f),
            (
                2,
                vec![
                    named("rsX", &x),
                    named("rsX", &[("GENE", "Y"), ("AA", "p.1")])
                ]
            )
        );
    }

    /// A later source whose every match at the record is lent to an earlier
    /// source's namesake lists no record of its own and still matched
    /// (`File.pm:339`, `:442`), so its statistics close over the match
    /// (`:346`) and are the short name's.
    #[test]
    fn a_source_lending_its_only_match_still_closes_its_statistics() {
        let dir = tempfile::tempdir().unwrap();
        let s = fixture(
            dir.path(),
            "s.vcf.gz",
            &HEADER,
            &[
                "21\t100\tCOSV1\tA\tC\t.\tPASS\tGENE=X",
                "21\t100\tCOSV2\tA\tC\t.\tPASS\tGENE=Z",
            ],
        );
        let b = fixture(
            dir.path(),
            "b.vcf.gz",
            &HEADER,
            &["21\t100\tCOSV1\tA\tC\t.\tPASS\tGENE=Y;AA=p.1"],
        );
        let cosmic = |file: &str| {
            format!(
                "file={file},short_name=COSMIC,format=vcf,type=exact,fields=GENE%AA,summary_stats=count"
            )
        };
        let sources = sources_of(&[&cosmic(&s), &cosmic(&b)], CustomRunOptions::default());
        let mut batch = vec![allele(100, 100, "A", "C")];
        annotate(&sources, &mut batch);
        let count = |n: f64| {
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(n))]
                    .into_iter()
                    .collect(),
            ))
        };
        assert_eq!(
            listed(&batch[0]),
            [
                named("COSV1", &[("GENE", "X"), ("AA", "p.1")]),
                named("COSV2", &[("GENE", "Z")]),
            ]
        );
        let [first, second] = batch[0].custom_annotations.as_slice() else {
            panic!("{:?}", batch[0].custom_annotations);
        };
        assert!(first.matched && first.stats == count(2.0), "{first:?}");
        assert!(
            second.records.is_empty() && second.matched && second.stats == count(1.0),
            "{second:?}"
        );
    }

    /// `VCF.pm:199-216`: a comma-separated value with one entry per alternate
    /// is split, each allele's record taking its own entry, and a value with
    /// one entry more drops its first, the reference's (`:203`); under
    /// `--custom_multi_allelic` the whole value stays on every allele. An
    /// `overlap` source rejoins the per-allele entries with `,` (`:255`).
    #[test]
    fn multi_allelic_keeps_the_whole_value_for_every_allele() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "multi.vcf.gz",
            &HEADER,
            &[
                "21\t100\ttwo\tA\tG,T\t.\tPASS\tAF=0.1,0.2",
                "21\t200\tthree\tA\tG,T\t.\tPASS\tAF=0.7,0.1,0.2",
                "21\t300\todd\tA\tG,T\t.\tPASS\tAF=0.1,0.2,0.3,0.4",
            ],
        );
        let values = |options: CustomRunOptions, kind: &str| -> Vec<Vec<String>> {
            let sources = sources(
                &format!("file={file},short_name=S,format=vcf,type={kind},fields=AF"),
                options,
            );
            let mut batch = Vec::new();
            for (record, pos) in [(1, 100), (2, 200), (3, 300)] {
                for alt in ["G", "T"] {
                    let mut variant = allele(pos, pos, "A", alt);
                    variant.input_record = record;
                    batch.push(variant);
                }
            }
            annotate(&sources, &mut batch);
            batch
                .iter()
                .map(|v| {
                    only(v)
                        .records
                        .iter()
                        .filter(|r| {
                            r.allele.is_none() || r.allele.as_deref() == Some(v.alt_allele())
                        })
                        .map(|r| text(&r.fields["AF"]).to_string())
                        .collect()
                })
                .collect()
        };
        let one = |s: &str| vec![s.to_string()];
        assert_eq!(
            values(CustomRunOptions::default(), "exact"),
            [
                one("0.1"),
                one("0.2"),
                one("0.1"),
                one("0.2"),
                one("0.1,0.2,0.3,0.4"),
                one("0.1,0.2,0.3,0.4"),
            ]
        );
        let multi_allelic = CustomRunOptions {
            multi_allelic: true,
            suppress_filter: false,
        };
        assert_eq!(
            values(multi_allelic, "exact"),
            [
                one("0.1,0.2"),
                one("0.1,0.2"),
                one("0.7,0.1,0.2"),
                one("0.7,0.1,0.2"),
                one("0.1,0.2,0.3,0.4"),
                one("0.1,0.2,0.3,0.4"),
            ]
        );
        assert_eq!(
            values(CustomRunOptions::default(), "overlap"),
            [
                one("0.1,0.2"),
                one("0.1,0.2"),
                one("0.1,0.2"),
                one("0.1,0.2"),
                one("0.1,0.2,0.3,0.4"),
                one("0.1,0.2,0.3,0.4"),
            ]
        );
    }

    /// `VCF.pm:196-199`: a file whose header says `##source=ClinVar` keeps
    /// every comma-separated value whole, as `--custom_multi_allelic` does.
    #[test]
    fn clinvar_source_header_never_splits() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "clinvar.vcf.gz",
            &[
                "##fileformat=VCFv4.1",
                "##source=ClinVar",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            ],
            &["21\t100\t12\tA\tG,T\t.\t.\tCLNSIG=Benign,Pathogenic;AF=0.1,0.2"],
        );
        let sources = sources(
            &format!("file={file},short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%AF"),
            CustomRunOptions::default(),
        );
        let mut batch = vec![allele(100, 100, "A", "G"), allele(100, 100, "A", "T")];
        annotate(&sources, &mut batch);
        for variant in &batch {
            let records = &only(variant).records;
            assert_eq!(records.len(), 2);
            for record in records {
                assert_eq!(text(&record.fields["CLNSIG"]), "Benign,Pathogenic");
                assert_eq!(text(&record.fields["AF"]), "0.1,0.2");
            }
        }
    }

    /// `VCF.pm:219-223`: a source with `fields` adds the record's FILTER
    /// entries, split at `;`, unless `--custom_suppress_filter` is on and
    /// `fields` does not name `FILTER`; a source without `fields` adds none.
    #[test]
    fn suppress_filter_drops_filter_unless_fields_name_it() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "filter.vcf.gz",
            &HEADER,
            &[
                "21\t100\tpass\tA\tG\t.\tPASS\tAF=0.1",
                "21\t200\tfail\tA\tG\t.\tq10;s50\tAF=0.2",
            ],
        );
        let filters = |fields: &str, suppress: bool| -> Vec<Option<CustomFieldValue>> {
            let sources = sources(
                &format!("file={file},short_name=S,format=vcf,type=exact{fields}"),
                CustomRunOptions {
                    multi_allelic: false,
                    suppress_filter: suppress,
                },
            );
            let mut batch = vec![allele(100, 100, "A", "G"), allele(200, 200, "A", "G")];
            batch[1].input_record = 2;
            annotate(&sources, &mut batch);
            batch
                .iter()
                .map(|v| only(v).records[0].fields.get("FILTER").cloned())
                .collect()
        };
        let list = |entries: &[&str]| {
            Some(CustomFieldValue::List(
                entries.iter().map(|e| e.to_string()).collect(),
            ))
        };
        assert_eq!(
            filters(",fields=AF", false),
            [list(&["PASS"]), list(&["q10", "s50"])]
        );
        assert_eq!(filters(",fields=AF", true), [None, None]);
        assert_eq!(
            filters(",fields=AF%FILTER", true),
            [list(&["PASS"]), list(&["q10", "s50"])]
        );
        assert_eq!(filters("", false), [None, None]);
    }

    /// `BaseVCF4.pm:86-163` on the header: a `##INFO` line's `Description` is
    /// kept by its `ID`, its quotes gone and its commas kept (`:139-150`), a
    /// later line with the same `ID` replacing the earlier, a line without a
    /// `Description` adding none; `##source=` is kept whole. The split at `=`
    /// keeps the piece after the first (`:146`), the angle brackets go
    /// wherever they stand (`:135`), and only the last quoted text keeps its
    /// commas. The output reads the descriptions for a file named `.vcf`,
    /// `.vcf.gz` or `.vcf.bgz` alone (`OutputFactory/BaseTab.pm:86`).
    #[test]
    fn info_header_reader_takes_descriptions_and_source() {
        let metadata = VcfMetadata::parse_lines([
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, per alternate\">",
            "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">",
            "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP member\">",
            "##INFO=<ID=NODESC,Number=1,Type=String>",
            "##INFO=<ID=EQ,Number=1,Type=String,Description=\"left=right\">",
            "##INFO=<ID=ANGLE,Number=1,Type=String,Description=\"a <b> c\">",
            "##INFO=<ID=TWO,Number=1,Type=String,Description=\"A, B\",Source=\"x, y\">",
            "##source=dbSNP",
            "##source=ClinVar",
            "## contig=<ID=21>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ])
        .unwrap();
        assert_eq!(metadata.source.as_deref(), Some("ClinVar"));
        let descriptions: Vec<(&str, &str)> = metadata
            .info_descriptions
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        assert_eq!(
            descriptions,
            [
                ("AF", "Allele frequency, per alternate"),
                ("DB", "dbSNP member"),
                ("EQ", "left"),
                ("ANGLE", "a b c"),
                ("TWO", "A"),
            ]
        );
        assert_eq!(VcfMetadata::parse_lines(["##source"]).unwrap().source, None);

        let dir = tempfile::tempdir().unwrap();
        let header = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ];
        let record = ["21\t100\trs1\tA\tG\t.\tPASS\tAF=0.1"];
        for (name, applies) in [
            ("sites.vcf.gz", true),
            ("sites.VCF.BGZ", true),
            ("sites.txt.gz", false),
        ] {
            let file = fixture(dir.path(), name, &header, &record);
            let sources = sources(
                &format!("file={file},format=vcf,type=exact,fields=AF"),
                CustomRunOptions::default(),
            );
            let source = sources.iter().next().unwrap();
            assert_eq!(
                source.field_descriptions().map(|d| d["AF"].as_str()),
                applies.then_some("Allele frequency"),
                "{name}"
            );
        }
    }

    /// The header is read as the bytes ensembl-io reads (`BaseVCF4.pm:86`): a
    /// byte that is not UTF-8 (Latin-1 `é`, 0xE9) in a description leaves the
    /// source's header readable, the byte written as U+FFFD in the
    /// description Ensembl VEP prints as it stands.
    #[test]
    fn a_header_that_is_not_utf8_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_bgzf_bytes(
            dir.path(),
            "sites.vcf.gz",
            b"##fileformat=VCFv4.2\n##source=caf\xE9\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Fr\xE9quence\">\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\trs1\tA\tG\t.\tPASS\tAF=0.1\n",
        );
        let metadata = VcfMetadata::read(&file).unwrap();
        assert_eq!(metadata.source.as_deref(), Some("caf\u{FFFD}"));
        assert_eq!(
            metadata.info_descriptions.get("AF").map(String::as_str),
            Some("Fr\u{FFFD}quence")
        );
    }

    /// Ensembl VEP opens a `--custom` VCF through the parser of its input
    /// (`VCF4Tabix.pm:55-58` reads every header line with `read_metadata`), so
    /// a `##fileformat` line naming a major version other than 4 ends the run
    /// as the source opens, with `BaseVCF4.pm:125`'s text, and a value with no
    /// `<digits>.<digits>` with `:129`'s; a `VCFv4.5` line and a file with no
    /// `##fileformat` line are read (measured on 116.2 with each shape as the
    /// first header line of the sites file: exit 2 and no output for `VCFv3.3`
    /// and `VCFv4`, the annotated row for `VCFv4.5` and the line's absence).
    /// htslib hands each header line over less its `\n` or `\r\n`, so
    /// `##fileformat=\r\n` is read as an empty value and `##fileformat=\r\r\n`
    /// refused, its value the second `\r` (both measured on 116.2).
    #[test]
    fn a_fileformat_line_outside_vcf_4_refuses_the_source_at_open() {
        let dir = tempfile::tempdir().unwrap();
        let info = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let record = ["21\t100\trs1\tA\tG\t.\tPASS\tAF=0.1"];
        let open = |name: &str, header: &[&str]| {
            let file = fixture(dir.path(), name, header, &record);
            let spec =
                parse_custom(&format!("file={file},format=vcf,type=exact,fields=AF")).unwrap();
            (
                file,
                VcfSource::open(spec, Arc::default(), CustomRunOptions::default()).map(|_| ()),
            )
        };

        let (file, refused) = open("v33.vcf.gz", &["##fileformat=VCFv3.3", info, columns]);
        assert_eq!(
            refused,
            Err(CustomError::Open {
                file: file.clone(),
                reason: "##fileformat=VCFv3.3: The VCF file format version 3.3 is not \
                         compatible with the parser version (VCF v4.2)"
                    .to_string(),
            })
        );
        assert_eq!(
            refused.unwrap_err().to_string(),
            format!(
                "cannot open the --custom source {file}: ##fileformat=VCFv3.3: The VCF file \
                 format version 3.3 is not compatible with the parser version (VCF v4.2)"
            )
        );
        let (file, refused) = open("v4.vcf.gz", &["##fileformat=VCFv4", info, columns]);
        assert_eq!(
            refused,
            Err(CustomError::Open {
                file,
                reason: "##fileformat=VCFv4: The script can't read the VCF file format \
                         version of 'fileformat'"
                    .to_string(),
            })
        );
        // The check runs on every header line, not the first alone (`:55-58`).
        let (file, refused) = open("late.vcf.gz", &[info, "##fileformat=VCFv5.0", columns]);
        assert_eq!(
            refused,
            Err(CustomError::Open {
                file,
                reason: "##fileformat=VCFv5.0: The VCF file format version 5.0 is not \
                         compatible with the parser version (VCF v4.2)"
                    .to_string(),
            })
        );

        let (file, refused) = open("cr.vcf.gz", &["##fileformat=\r\r", info, columns]);
        assert_eq!(
            refused,
            Err(CustomError::Open {
                file,
                reason: "##fileformat=\r: The script can't read the VCF file format \
                         version of 'fileformat'"
                    .to_string(),
            })
        );

        assert_eq!(open("v42.vcf.gz", &[HEADER[0], info, columns]).1, Ok(()));
        assert_eq!(
            open("v45.vcf.gz", &["##fileformat=VCFv4.5", info, columns]).1,
            Ok(())
        );
        assert_eq!(open("none.vcf.gz", &[info, columns]).1, Ok(()));
        assert_eq!(
            open(
                "crlf.vcf.gz",
                &["##fileformat=VCFv4.2\r", "##fileformat=\r", columns]
            )
            .1,
            Ok(())
        );
    }

    /// The record `T TCT CLNSIG=Pathogenic` and `T TTC CLNSIG=Benign` at
    /// 21:25593910 against the input `T TCT`, which inserts `CT` into a `CT`
    /// repeat, and `T TCT,TTC`: each alternate keeps the record whose
    /// alleles it matches, keyed on the allele as the input gives it, so a
    /// 3' shift of the printed allele on a coding transcript (`CT` printed
    /// `TC`) changes nothing here. Expected values: release/116.3
    /// `t/OutputFactory.t:1894-1951` (`c511af47`), `ClinVar_CLNSIG`
    /// `Pathogenic` under shift 0 and 1, and `Pathogenic` for the first
    /// alternate and `Benign` for the second of the multi-allelic record.
    #[test]
    fn shifted_coding_rows_keep_their_own_allele_records() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "test_shift_3prime.vcf.gz",
            &[
                "##fileformat=VCFv4.1",
                "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            ],
            &[
                "21\t25593910\tins_ct\tT\tTCT\t.\tPASS\tCLNSIG=Pathogenic",
                "21\t25593910\tins_tc\tT\tTTC\t.\tPASS\tCLNSIG=Benign",
            ],
        );
        let sources = sources(
            &format!("file={file},short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG"),
            CustomRunOptions::default(),
        );
        // `21 25593910 T TCT` chopped: the insertion `-/CT` at 25593911-25593910.
        let mut batch = vec![allele(25593911, 25593910, "-", "CT")];
        annotate(&sources, &mut batch);
        let records = &only(&batch[0]).records;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "ins_ct");
        assert_eq!(records[0].allele.as_deref(), Some(&b"CT"[..]));
        assert_eq!(text(&records[0].fields["CLNSIG"]), "Pathogenic");

        let mut batch = vec![
            allele(25593911, 25593910, "-", "CT"),
            allele(25593911, 25593910, "-", "TC"),
        ];
        annotate(&sources, &mut batch);
        let by_allele: Vec<(&[u8], &str)> = only(&batch[0])
            .records
            .iter()
            .map(|r| (r.allele.as_deref().unwrap(), text(&r.fields["CLNSIG"])))
            .collect();
        assert_eq!(
            by_allele,
            [(&b"CT"[..], "Pathogenic"), (&b"TC"[..], "Benign")]
        );
        assert_eq!(batch[0].custom_annotations, batch[1].custom_annotations);
    }

    /// `InputBuffer.pm:638-710` splits a multi-allelic record under `--minimal`
    /// into one feature per allele once the trimming moved one, and the rejoin
    /// (`OutputFactory.pm:2282-2340`) copies the later features' transcript
    /// variations alone, so the record's custom annotations are the first
    /// feature's: the record `ATG/ACG/A` at 100, split to `T/C` at 101 and
    /// `TG/-` at 101-102, lists the record its first allele matches and not
    /// the deletion its second would, and every allele carries that list. The
    /// same alleles left whole list both records, each keyed by the allele as
    /// the input gives it (`CG`, not the trimmed `C`). A record Ensembl VEP
    /// skips (`File.pm:318`, `vep_skip`) or one past `--max_sv_size` passes
    /// through unannotated, and a reference row shares its record's list
    /// without being matched as an allele.
    #[test]
    fn minimal_split_records_match_through_their_first_allele_alone() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture(
            dir.path(),
            "split.vcf.gz",
            &HEADER,
            &[
                "21\t100\tdel\tATG\tA\t.\tPASS\tDP=1",
                "21\t101\tsnv\tT\tC\t.\tPASS\tDP=2",
            ],
        );
        let sources = sources(
            &format!("file={file},short_name=S,format=vcf,type=exact,fields=DP"),
            CustomRunOptions::default(),
        );
        let listed = |variant: &InputVariant| -> Vec<(String, Vec<u8>)> {
            only(variant)
                .records
                .iter()
                .map(|r| (r.name.clone(), r.allele.clone().unwrap()))
                .collect()
        };

        let mut split = vec![allele(101, 101, "T", "C"), allele(101, 102, "TG", "-")];
        for variant in &mut split {
            variant.original_start = Some(101);
            variant.original_end = Some(102);
        }
        annotate(&sources, &mut split);
        assert_eq!(listed(&split[0]), [("snv".to_string(), b"C".to_vec())]);
        assert_eq!(split[0].custom_annotations, split[1].custom_annotations);

        let mut whole = vec![allele(101, 102, "TG", "CG"), allele(101, 102, "TG", "-")];
        annotate(&sources, &mut whole);
        assert_eq!(
            listed(&whole[0]),
            [
                ("del".to_string(), b"-".to_vec()),
                ("snv".to_string(), b"CG".to_vec()),
            ]
        );
        assert_eq!(whole[0].custom_annotations, whole[1].custom_annotations);

        let mut skipped = allele(101, 101, "T", "C");
        skipped.vep_skip = true;
        let mut oversize = allele(101, 101, "T", "C");
        oversize.oversize_sv = true;
        let mut batch = records(vec![skipped, oversize]);
        annotate(&sources, &mut batch);
        assert!(batch.iter().all(|v| v.custom_annotations.is_empty()));

        let alternate = allele(101, 101, "T", "C");
        let mut batch = vec![alternate.reference_allele_row(), alternate];
        annotate(&sources, &mut batch);
        assert_eq!(listed(&batch[1]), [("snv".to_string(), b"C".to_vec())]);
        assert_eq!(batch[0].custom_annotations, batch[1].custom_annotations);
    }
}
