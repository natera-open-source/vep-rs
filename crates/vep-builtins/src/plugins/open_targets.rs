// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! OpenTargets plugin: Open Targets Platform GWAS and QTL associations.
//!
//! Reads the Platform's tabix-indexed VEP export, a TSV whose first four
//! columns name the variant (`chromosome`, `position`, `referenceAllele`,
//! `alternateAllele`) and whose remaining columns describe one association per
//! line. Every selected column is written as `OpenTargets_<column>`, its values
//! `,`-joined across the records the variant matches, `NA` where a record has
//! none; the JSON format gets one object per record under `OpenTargets`
//! instead ([`InputVariant::plugin_records`]).
//!
//! Parameters are `key=value`: `file` (the export), `cols` (`all`, or a
//! `:`-separated list of header columns; the defaults are the GWAS columns
//! `gwasLocusToGeneScore`, `gwasGeneId` and `gwasDiseases` and the QTL columns
//! `qtlGeneId` and `qtlBiosampleName`, each family when its study types are
//! selected and its columns are in the header), `study_types` (`all`, `QTL`
//! for every study type ending in `qtl`, `gwas`, or a `:`-separated list of
//! such names) and `lead_variants_only` (`1`, the default, keeps the records
//! whose `isLead` is `true`).
//!
//! `OpenTargets.pm` is the plugin's counterpart in Ensembl/VEP_plugins
//! release/116; the Platform publishes GRCh38 coordinates only.

use std::borrow::Cow;
use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::allele_match::{self, KeyPosition};
use crate::tabix::{
    split_field_count, table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::variant::{InputVariant, PluginRecord};

const PLUGIN_NAME: &str = "OpenTargets";
const KEY_PREFIX: &str = "OpenTargets_";
const PARAMETERS: [&str; 4] = ["file", "cols", "study_types", "lead_variants_only"];
/// The export's leading columns, which identify the variant and are no output.
const VARIANT_COLUMNS: usize = 4;
const POS_COLUMN: usize = 1;
const REF_COLUMN: usize = 2;
const ALT_COLUMN: usize = 3;
const DEFAULT_GWAS_COLUMNS: [&str; 3] = ["gwasLocusToGeneScore", "gwasGeneId", "gwasDiseases"];
const DEFAULT_QTL_COLUMNS: [&str; 2] = ["qtlGeneId", "qtlBiosampleName"];
const IS_LEAD_KEY: &str = "OpenTargets_isLead";
const STUDY_TYPE_KEY: &str = "OpenTargets_studyType";
const L2G_KEY: &str = "OpenTargets_gwasLocusToGeneScore";
/// The text formats' value for a record without one (`_join_results`,
/// OpenTargets.pm:339).
const MISSING: &str = "NA";
/// A record is read from two bases before the variant's start to its end, so
/// a VCF-style record whose position is its anchor base is found for an
/// insertion or deletion (`run`, OpenTargets.pm:350).
const QUERY_LEAD: u64 = 2;

/// The region a variant's records are read from, as the prefetch queries it
/// and the annotation looks it up: [`QUERY_LEAD`] bases before the start, to
/// the end.
fn query_region(variant: &InputVariant, assembly: Option<Assembly>) -> (&str, u64, u64) {
    (
        table_contig(variant, assembly),
        variant.start.saturating_sub(QUERY_LEAD),
        variant.end,
    )
}

/// Perl's `split` without a limit, which reads the `cols` and `study_types`
/// lists: the pieces of `text` at `separator`, the trailing empty ones dropped
/// and the empty ones inside the list kept.
fn perl_split(text: &str, separator: char) -> Vec<&str> {
    let mut pieces: Vec<&str> = text.split(separator).collect();
    while pieces.last().is_some_and(|piece| piece.is_empty()) {
        pieces.pop();
    }
    pieces
}

/// The study types of one family (GWAS, or the QTL types) that pass the
/// `study_types` filter (`_get_studytypes`, OpenTargets.pm:126-209).
#[derive(Debug, Clone, PartialEq, Eq)]
enum TypeSet {
    /// None of the family's types.
    None,
    /// Every type of the family: `gwas`, or every name ending in `qtl`.
    Every,
    /// The named types.
    Named(Vec<String>),
}

impl TypeSet {
    fn is_empty(&self) -> bool {
        match self {
            TypeSet::None => true,
            TypeSet::Every => false,
            TypeSet::Named(names) => names.is_empty(),
        }
    }
}

/// The `study_types` selection: which GWAS and which QTL study types a record
/// may carry.
#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyTypes {
    gwas: TypeSet,
    qtl: TypeSet,
    /// Whether records are filtered on their `studyType` at all: Ensembl VEP
    /// filters when `study_types` was given and is not exactly `all`
    /// (OpenTargets.pm:359).
    filters: bool,
}

impl StudyTypes {
    /// Parses the `study_types` request against the families the header carries. An absent
    /// request, or one Perl reads as false (`0` or empty, `if($studytype_request_str)`,
    /// OpenTargets.pm:147), selects every family with columns, where Ensembl VEP selects
    /// both whatever the header carries (OpenTargets.pm:202-205); `all` selects every family
    /// with columns; `QTL` every QTL type; `gwas` and names ending in `qtl` are listed, and
    /// a listed QTL type beside `QTL` is a conflict, as is a family without columns.
    fn parse(
        request: Option<&str>,
        gwas_columns: bool,
        qtl_columns: bool,
        file: &str,
    ) -> Result<Self, PluginError> {
        let Some(request) = request.filter(|r| perl_true(r)) else {
            return Ok(Self {
                gwas: if gwas_columns {
                    TypeSet::Every
                } else {
                    TypeSet::None
                },
                qtl: if qtl_columns {
                    TypeSet::Every
                } else {
                    TypeSet::None
                },
                filters: false,
            });
        };
        let no_columns = |family: &str, token: &str| {
            PluginError::Init(format!(
                "{PLUGIN_NAME}: study type '{token}' in '{request}' matches {family} study \
                 types but the header of {file} has no {family} columns"
            ))
        };
        let mut gwas = TypeSet::None;
        let mut qtl = TypeSet::None;
        for token in perl_split(request, ':') {
            match token {
                "all" => {
                    if gwas_columns {
                        gwas = TypeSet::Every;
                    }
                    if qtl_columns {
                        qtl = TypeSet::Every;
                    }
                    break;
                }
                "QTL" => {
                    if matches!(&qtl, TypeSet::Named(names) if !names.is_empty()) {
                        return Err(PluginError::Init(format!(
                            "{PLUGIN_NAME}: study type 'QTL' in '{request}' conflicts with the \
                             QTL study types listed before it"
                        )));
                    }
                    if !qtl_columns {
                        return Err(no_columns("QTL", token));
                    }
                    qtl = TypeSet::Every;
                }
                "gwas" => {
                    if !gwas_columns {
                        return Err(no_columns("GWAS", token));
                    }
                    if let TypeSet::Named(names) = &mut gwas {
                        names.push(token.to_string());
                    } else {
                        gwas = TypeSet::Named(vec![token.to_string()]);
                    }
                }
                _ if token.ends_with("qtl") => match &mut qtl {
                    TypeSet::Every => {
                        return Err(PluginError::Init(format!(
                            "{PLUGIN_NAME}: study type '{token}' in '{request}' conflicts with \
                             the QTL before it, which already selects every QTL study type"
                        )));
                    }
                    TypeSet::Named(names) => names.push(token.to_string()),
                    TypeSet::None => {
                        if !qtl_columns {
                            return Err(no_columns("QTL", token));
                        }
                        qtl = TypeSet::Named(vec![token.to_string()]);
                    }
                },
                _ => {
                    return Err(PluginError::Init(format!(
                        "{PLUGIN_NAME}: unknown study type '{token}' in '{request}'; a study \
                         type is all, QTL, gwas or a name ending in qtl"
                    )));
                }
            }
        }
        Ok(Self {
            gwas,
            qtl,
            filters: request != "all",
        })
    }

    /// Whether a record of `study_type` passes the filter.
    fn accepts(&self, study_type: &str) -> bool {
        let in_family = |set: &TypeSet, every: bool| match set {
            TypeSet::None => false,
            TypeSet::Every => every,
            TypeSet::Named(names) => names.iter().any(|name| name == study_type),
        };
        in_family(&self.gwas, study_type == "gwas")
            || in_family(&self.qtl, study_type.ends_with("qtl"))
    }
}

/// The value of `record` under the plugin's key at `index` (the key's place
/// among the columns after the variant columns); `None` past the line's last
/// filled cell, `present` ([`split_field_count`]): a record whose last columns
/// are empty has no value there (`undef`, `null` in the JSON format), where an
/// empty cell before a filled one is the empty string.
fn cell(record: &TabixRecord, present: usize, index: usize) -> Option<&str> {
    let column = VARIANT_COLUMNS + index;
    if column < present {
        record.columns.get(column).map(String::as_str)
    } else {
        None
    }
}

/// A record's alleles and start as Ensembl VEP compares them (`parse_data`,
/// OpenTargets.pm:394-411): a VCF-style insertion or deletion, whose alleles
/// differ in length and share their first base, loses that base and starts one
/// base later, an emptied allele written `-`. The start is the position cell's
/// text, the key's text until the trim inside the match moves it
/// ([`KeyPosition::Cell`]); the base lost here steps it as `$start++` (`:405`)
/// steps a cell no arithmetic has read ([`allele_match::perl_increment`]:
/// `6000010.0` is `6000011`, `06000010` is `06000011`, `abc` is `abd`). `None`
/// for a record whose reference or alternate cell lies past the line's last
/// filled cell, the undefined value `split /\t/` leaves (`:396`,
/// [`allele_match::defined_alleles`]); an empty cell before a filled one is the
/// defined empty string, compared as `-`.
fn record_alleles(record: &TabixRecord) -> Option<(&[u8], &[u8], Cow<'_, str>)> {
    let (reference, alt) = allele_match::defined_alleles(record, REF_COLUMN, ALT_COLUMN)?;
    let (reference, alt) = (reference.as_bytes(), alt.as_bytes());
    let cell = record.columns.get(POS_COLUMN).map_or("", String::as_str);
    if reference.len() != alt.len() && reference.first().is_some_and(|b| alt.first() == Some(b)) {
        return Some((
            or_dash(&reference[1..]),
            or_dash(&alt[1..]),
            Cow::Owned(allele_match::perl_increment(cell)),
        ));
    }
    Some((reference, alt, Cow::Borrowed(cell)))
}

/// `-`, the empty allele's spelling, for an allele with no base left.
fn or_dash(allele: &[u8]) -> &[u8] {
    if allele.is_empty() {
        b"-"
    } else {
        allele
    }
}

/// The reference and alternate allele with the bases they share trimmed away
/// and the start moved past every base trimmed on the left, `-` for an allele
/// trimmed to nothing (`trim_sequences`, ensembl-variation
/// Utils/Sequence.pm:965-1038, on the forward strand): the left end first, or
/// the right end first.
fn trimmed<'a>(
    reference: &[u8],
    alt: &[u8],
    start: KeyPosition<'a>,
    right_first: bool,
) -> (Vec<u8>, Vec<u8>, KeyPosition<'a>) {
    let mut reference = reference;
    let mut alt = alt;
    let mut start = start;
    let trim_left = |reference: &mut &[u8], alt: &mut &[u8], start: &mut KeyPosition<'a>| {
        while let (Some(r), Some(a)) = (reference.first(), alt.first()) {
            if r != a {
                break;
            }
            *reference = &reference[1..];
            *alt = &alt[1..];
            *start = start.moved();
        }
    };
    let trim_right = |reference: &mut &[u8], alt: &mut &[u8]| {
        while let (Some(r), Some(a)) = (reference.last(), alt.last()) {
            if r != a {
                break;
            }
            *reference = &reference[..reference.len() - 1];
            *alt = &alt[..alt.len() - 1];
        }
    };
    if right_first {
        trim_right(&mut reference, &mut alt);
        trim_left(&mut reference, &mut alt, &mut start);
    } else {
        trim_left(&mut reference, &mut alt, &mut start);
        trim_right(&mut reference, &mut alt);
    }
    (or_dash(reference).to_vec(), or_dash(alt).to_vec(), start)
}

/// The trimmed forms of a pair of alleles: left end first, and right end first
/// when either allele is longer than one base (`_get_trim_directions`,
/// Utils/Sequence.pm:1218-1225).
fn trimmed_forms<'a>(
    reference: &[u8],
    alt: &[u8],
    start: KeyPosition<'a>,
) -> Vec<(Vec<u8>, Vec<u8>, KeyPosition<'a>)> {
    let mut forms = vec![trimmed(reference, alt, start, false)];
    if reference.len() > 1 || alt.len() > 1 {
        forms.push(trimmed(reference, alt, start, true));
    }
    forms
}

/// The alternate alleles of the whole input record on the forward strand: the
/// plugin matches a data record against every alternate allele of the record
/// (`alt_alleles`, OpenTargets.pm:349), so a record matched by one allele of a
/// multi-allelic record is written on every allele's rows.
fn record_alternate_alleles(variant: &InputVariant) -> Vec<Vec<u8>> {
    // The borrowed form is the one sign nothing was complemented: a `-` or a
    // palindrome such as `AT` reads the same on both strands.
    let (forward_ref, _) = variant.forward_strand_alleles();
    let complemented = matches!(forward_ref, Cow::Owned(_));
    variant
        .record_allele_string()
        .split('/')
        .skip(1)
        .map(|allele| {
            if complemented {
                vep_core::codon::reverse_complement(allele.as_bytes())
            } else {
                allele.as_bytes().to_vec()
            }
        })
        .collect()
}

/// Whether a record's alleles name an allele of the input record
/// (`get_matched_variant_alleles`, Utils/Sequence.pm:1144-1216): a trimmed
/// form of the record's alleles and start equals one of the variant's, the
/// variant's alleles read on the forward strand and the starts compared as
/// the key holds them ([`KeyPosition::equals`]). A record whose position is
/// Perl-false as `parse_data` leaves it matches nothing (`:1173`), whatever its
/// allele cells hold. `None` for a record with no reference or alternate allele
/// cell ([`record_alleles`]) and a Perl-true position, which the match dies on
/// before any comparison (`trim_sequences`, `:968`, reached from `:1230`).
fn alleles_match(
    variant: &InputVariant,
    alternates: &[Vec<u8>],
    record: &TabixRecord,
) -> Option<bool> {
    // `!$b->{pos}` (Sequence.pm:1173) returns the match before `trim_sequences`
    // dies on an undefined allele (`:1230`), and reads the cell after
    // `parse_data`'s move (OpenTargets.pm:405), which steps `0` and the empty
    // cell to 1: an undefined allele never moves the cell, so a Perl-false cell
    // beside one ends the match.
    let position_cell = record.columns.get(POS_COLUMN).map_or("", String::as_str);
    let Some((reference, alt, position)) = record_alleles(record) else {
        return if perl_true(position_cell) {
            None
        } else {
            Some(false)
        };
    };
    if !perl_true(&position) {
        return Some(false);
    }
    let start = KeyPosition::Cell {
        text: &position,
        moved: 0,
    };
    let (variant_ref, _) = variant.forward_strand_alleles();
    let record_forms = trimmed_forms(reference, alt, start);
    Some(alternates.iter().any(|variant_alt| {
        let variant_forms = trimmed_forms(
            &variant_ref,
            variant_alt,
            KeyPosition::Number(variant.start),
        );
        record_forms.iter().any(|(reference, alt, start)| {
            variant_forms
                .iter()
                .any(|(variant_reference, variant_alt, variant_start)| {
                    variant_reference == reference
                        && variant_alt == alt
                        && start.equals(*variant_start)
                })
        })
    }))
}

/// Open Targets Platform association plugin.
pub struct OpenTargetsPlugin {
    /// The export, opened by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    /// The export's path as given, which a failure on one of its rows names.
    file: String,
    /// The file's name, which every field description cites.
    file_name: String,
    /// Every header column after the variant columns, prefixed, in file order:
    /// the keys of a record.
    record_keys: Vec<String>,
    /// The keys written to the output, in the order `cols` gave them.
    selected: Vec<String>,
    /// Where `isLead` and `studyType` sit among `record_keys`.
    is_lead: Option<usize>,
    study_type: Option<usize>,
    lead_variants_only: bool,
    study_types: StudyTypes,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for OpenTargetsPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl OpenTargetsPlugin {
    /// Create an uninitialized OpenTargets plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            file: String::new(),
            file_name: String::new(),
            record_keys: Vec::new(),
            selected: Vec::new(),
            is_lead: None,
            study_type: None,
            lead_variants_only: true,
            study_types: StudyTypes {
                gwas: TypeSet::None,
                qtl: TypeSet::None,
                filters: false,
            },
            assembly: None,
        }
    }

    /// The output key of a header column: `OpenTargets_` before it, unless
    /// the column already carries the prefix (`_prefix_cols`,
    /// OpenTargets.pm:98-103).
    fn key_of(column: &str) -> String {
        if column.starts_with(KEY_PREFIX) {
            column.to_string()
        } else {
            format!("{KEY_PREFIX}{column}")
        }
    }

    /// The columns to write: `cols` as given, a name listed twice written once
    /// (`:308` and `:318` key the header and the row by a hash over the list),
    /// or, when it is absent or Perl reads it as false (`0` or empty,
    /// `if(!$cols)`, OpenTargets.pm:271), the defaults of each study family that
    /// is selected and has columns (`:269-282`).
    fn columns_to_write(
        cols: Option<&str>,
        columns: &[String],
        study_types: &StudyTypes,
        gwas_columns: bool,
        qtl_columns: bool,
        file: &str,
    ) -> Result<Vec<String>, PluginError> {
        let wanted: Vec<String> = match cols.filter(|c| perl_true(c)) {
            Some("all") => return Ok(columns.to_vec()),
            Some(cols) => perl_split(cols, ':')
                .into_iter()
                .map(str::to_string)
                .collect(),
            None => {
                let mut names = Vec::new();
                if !study_types.gwas.is_empty() && gwas_columns {
                    names.extend(DEFAULT_GWAS_COLUMNS.map(str::to_string));
                }
                if !study_types.qtl.is_empty() && qtl_columns {
                    names.extend(DEFAULT_QTL_COLUMNS.map(str::to_string));
                }
                names
            }
        };
        let missing: Vec<&str> = wanted
            .iter()
            .filter(|name| !columns.contains(name))
            .map(|name| if name.is_empty() { "''" } else { name.as_str() })
            .collect();
        if !missing.is_empty() {
            return Err(PluginError::Init(format!(
                "{PLUGIN_NAME}: columns not found in the header of {file}: {}",
                missing.join(", ")
            )));
        }
        let mut once: Vec<String> = Vec::with_capacity(wanted.len());
        for name in wanted {
            if !once.contains(&name) {
                once.push(name);
            }
        }
        Ok(once)
    }

    /// Whether `record` passes the `lead_variants_only` and `study_types` filters
    /// (OpenTargets.pm:352-372).
    fn keeps(&self, record: &TabixRecord) -> bool {
        let present = split_field_count(&record.columns);
        if self.lead_variants_only
            && self
                .is_lead
                .is_none_or(|index| cell(record, present, index) != Some("true"))
        {
            return false;
        }
        if self.study_types.filters {
            return self.study_type.is_some_and(|index| {
                cell(record, present, index)
                    .is_some_and(|study_type| self.study_types.accepts(study_type))
            });
        }
        true
    }

    /// The records of `candidates` that annotate `variant`, in file order, as
    /// `run` selects them (OpenTargets.pm:345-392): the `lead_variants_only`
    /// and `study_types` filters drop their records first (`:352-372`,
    /// [`Self::keeps`]), then every record left is put to the allele match
    /// (`:375-390`), which dies on a record with no reference or alternate
    /// allele cell before any comparison (`trim_sequences`,
    /// Utils/Sequence.pm:968). Such a record fails the variant wherever it
    /// lies among the records that match, while one the filters drop never
    /// reaches the match.
    fn matched<'a>(
        &self,
        variant: &InputVariant,
        candidates: Vec<&'a TabixRecord>,
    ) -> Result<Vec<&'a TabixRecord>, PluginError> {
        let alternates = record_alternate_alleles(variant);
        let mut matched = Vec::new();
        for record in candidates.into_iter().filter(|record| self.keeps(record)) {
            match alleles_match(variant, &alternates, record) {
                Some(true) => matched.push(record),
                Some(false) => {}
                None => {
                    return Err(allele_match::undefined_allele_row(
                        PLUGIN_NAME,
                        &self.file,
                        record,
                    ))
                }
            }
        }
        Ok(matched)
    }

    /// Writes the matched records onto `variant`: each selected key's values
    /// `,`-joined in file order for the text formats, `NA` for a record without
    /// one, and one record object per match for the JSON format
    /// (`_join_results`, OpenTargets.pm:322-343).
    fn write(&self, variant: &mut InputVariant, matched: &[&TabixRecord]) {
        let selected: Vec<usize> = self
            .selected
            .iter()
            .filter_map(|key| self.record_keys.iter().position(|k| k == key))
            .collect();
        let mut joined: Vec<String> = vec![String::new(); selected.len()];
        let mut records: Vec<PluginRecord> = Vec::with_capacity(matched.len());
        for (n, record) in matched.iter().enumerate() {
            let present = split_field_count(&record.columns);
            let mut object = PluginRecord::with_capacity(selected.len());
            for (slot, &index) in selected.iter().enumerate() {
                let value = cell(record, present, index);
                if n > 0 {
                    joined[slot].push(',');
                }
                joined[slot].push_str(value.filter(|v| !v.is_empty()).unwrap_or(MISSING));
                object.insert(self.record_keys[index].clone(), value.map(str::to_string));
            }
            records.push(object);
        }
        for (slot, &index) in selected.iter().enumerate() {
            variant.plugin_data.insert(
                self.record_keys[index].clone(),
                std::mem::take(&mut joined[slot]),
            );
        }
        variant
            .plugin_records
            .insert(PLUGIN_NAME.to_string(), records);
    }
}

impl BuiltinPlugin for OpenTargetsPlugin {
    fn name(&self) -> &str {
        PLUGIN_NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`OpenTargets.pm:291-292`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        &RowKind::ALL
    }

    /// Every selected key is "column from <file>"; the locus-to-gene score
    /// says what it is first (`get_header_info`, OpenTargets.pm:295-311).
    fn header_info(&self) -> Vec<(String, String)> {
        let column = format!("column from {}", self.file_name);
        self.selected
            .iter()
            .map(|key| {
                let description = if key == L2G_KEY {
                    format!(
                        "Locus-to-gene (L2G) scores to predict causal genes at GWAS loci; {column}"
                    )
                } else {
                    column.clone()
                };
                (key.clone(), description)
            })
            .collect()
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(PLUGIN_NAME, &PARAMETERS, params)?;
        let file = params
            .get("file")
            .filter(|f| !f.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{PLUGIN_NAME}: no file specified; give file=<path to the Open Targets \
                     export>"
                ))
            })?;
        self.lead_variants_only = match params.get("lead_variants_only") {
            None | Some("1") => true,
            Some("" | "0") => false,
            Some(other) => {
                return Err(PluginError::Init(format!(
                    "{PLUGIN_NAME}: lead_variants_only must be 0 or 1, not '{other}'"
                )))
            }
        };

        let config = TabixConfig {
            file_path: file.into(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };
        // A row keeps the whitespace before its line end, which Ensembl VEP's
        // `split /\t/` leaves on the line's last cell (`parse_data`,
        // OpenTargets.pm:396), so a whitespace-only last cell is a value.
        let store = open_best_store_with(
            Path::new(file),
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )
        .map_err(|e| match e {
            PluginError::Init(message) => PluginError::Init(format!("{PLUGIN_NAME}: {message}")),
            other => other,
        })?;
        let header = store.header().ok_or_else(|| {
            PluginError::Init(format!("{PLUGIN_NAME}: no header line found in {file}"))
        })?;
        let columns: Vec<String> = header.iter().skip(VARIANT_COLUMNS).cloned().collect();
        if columns.is_empty() {
            return Err(PluginError::Init(format!(
                "{PLUGIN_NAME}: the header of {file} has no columns after the four variant \
                 columns"
            )));
        }
        let gwas_columns = columns.iter().any(|c| c.starts_with("gwas"));
        let qtl_columns = columns.iter().any(|c| c.starts_with("qtl"));
        let study_types =
            StudyTypes::parse(params.get("study_types"), gwas_columns, qtl_columns, file)?;
        if study_types.gwas.is_empty() && study_types.qtl.is_empty() {
            return Err(PluginError::Init(format!(
                "{PLUGIN_NAME}: {file} carries neither GWAS nor QTL data: its header has no \
                 gwas or qtl columns"
            )));
        }
        let selected = Self::columns_to_write(
            params.get("cols"),
            &columns,
            &study_types,
            gwas_columns,
            qtl_columns,
            file,
        )?;

        self.record_keys = columns.iter().map(|c| Self::key_of(c)).collect();
        self.selected = selected.iter().map(|c| Self::key_of(c)).collect();
        self.is_lead = self.record_keys.iter().position(|k| k == IS_LEAD_KEY);
        self.study_type = self.record_keys.iter().position(|k| k == STUDY_TYPE_KEY);
        if self.lead_variants_only && self.is_lead.is_none() {
            return Err(PluginError::Init(format!(
                "{PLUGIN_NAME}: lead_variants_only=1 needs an isLead column, which the header \
                 of {file} lacks; pass lead_variants_only=0"
            )));
        }
        if study_types.filters && self.study_type.is_none() {
            return Err(PluginError::Init(format!(
                "{PLUGIN_NAME}: study_types needs a studyType column, which the header of \
                 {file} lacks"
            )));
        }
        self.study_types = study_types;
        self.file = file.to_string();
        self.file_name = Path::new(file)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| file.to_string());

        debug!(
            file,
            selected = ?self.selected,
            lead_variants_only = self.lead_variants_only,
            "OpenTargets plugin initialized"
        );
        self.store = Some(store);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(None);
        }
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| self.annotates_variant(v))
            .map(|v| query_region(v, self.assembly))
            .collect();
        Ok(Some(PrefetchData::new(store.query_batch(&regions)?)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;
        if batch.is_empty() {
            return Ok(());
        }
        let mut failure = None;
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let (chr, start, end) = query_region(variant, self.assembly);
            let candidates = batch.lookup(chr, start, end);
            if candidates.is_empty() {
                continue;
            }
            match self.matched(variant, candidates) {
                Ok(matched) if !matched.is_empty() => self.write(variant, &matched),
                Ok(_) => {}
                Err(error) => {
                    failure.get_or_insert(error);
                }
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `annotate`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use vep_core::coordinate::Strand;

    /// The export's header: the four variant columns, then the association
    /// columns of the Platform's VEP export.
    const HEADER: &str = "chromosome\tposition\treferenceAllele\talternateAllele\tvariantId\
                          \tstudyId\tstudyType\tstudyLocusId\tpValueMantissa\tpValueExponent\
                          \tbeta\tis95CredibleSet\tis99CredibleSet\tposteriorProbability\
                          \tfinemappingMethod\tisLead\tgwasGeneId\tgwasLocusToGeneScore\
                          \tgwasDiseases\tqtlGeneId\tqtlBiosampleName";

    /// Records shaped like the export's: two GWAS associations of one variant
    /// (one a lead variant, one not) whose QTL cells trail empty, an eQTL and a
    /// pQTL association, a GWAS record of a VCF-style deletion and one of a
    /// VCF-style insertion.
    const RECORDS: [&str; 6] = [
        "21\t5030000\tA\tG\t21_5030000_A_G\tGCST1\tgwas\tloc1\t6.331\t-11\t-0.448027\ttrue\ttrue\
         \t1.0\tPICS\ttrue\tENSG00000000001\t0.5500299334526062\tEFO_0004527\t\t",
        "21\t5030000\tA\tG\t21_5030000_A_G\tGCST2\tgwas\tloc2\t1.275\t-10\t-0.437719\ttrue\ttrue\
         \t1.0\tPICS\tfalse\tENSG00000000002\t0.61\tEFO_0000001\t\t",
        "21\t5030000\tA\tG\t21_5030000_A_G\tgtex_ge_cerebellum\teqtl\tloc3\t6.103\t-7\t1.05262\
         \ttrue\ttrue\t0.95\tSuSie\ttrue\t\t\t\tENSG00000000003\tcerebellum",
        "21\t5030001\tC\tT\t21_5030001_C_T\tukb_ppp_liver\tpqtl\tloc4\t2.5\t-9\t0.3\ttrue\ttrue\
         \t0.9\tSuSie\ttrue\t\t\t\tENSG00000000004\tliver",
        "21\t5030010\tAT\tA\t21_5030010_AT_A\tGCST3\tgwas\tloc5\t3.0\t-8\t0.1\ttrue\ttrue\t1.0\
         \tPICS\ttrue\tENSG00000000005\t0.7\tEFO_0000002\t\t",
        "21\t5030020\tG\tGCA\t21_5030020_G_GCA\tGCST4\tgwas\tloc6\t4.0\t-12\t0.2\ttrue\ttrue\
         \t1.0\tPICS\ttrue\tENSG00000000006\t0.8\tEFO_0000003\t\t",
    ];

    /// The fixture as the plugin's instructions prepare the export: the header
    /// commented with `#`, indexed on the chromosome and position columns.
    fn export(dir: &Path) -> String {
        let header = format!("#{HEADER}");
        let mut lines = vec![header.as_str()];
        lines.extend(RECORDS);
        write_tabix_fixture(
            dir,
            "open_targets_vep.tsv.gz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string()
    }

    /// The export as `name`, with `rows` spliced in before the record at
    /// `index` of [`RECORDS`] (after the last one at `RECORDS.len()`).
    fn export_with(dir: &Path, name: &str, index: usize, rows: &[&str]) -> String {
        let header = format!("#{HEADER}");
        let mut lines = vec![header.as_str()];
        lines.extend(&RECORDS[..index]);
        lines.extend(rows);
        lines.extend(&RECORDS[index..]);
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    /// `file=<file>` followed by `params`, the argument list `init` takes.
    fn init_params(file: &str, params: &[&str]) -> Vec<String> {
        let mut all = vec![format!("file={file}")];
        all.extend(params.iter().map(|p| p.to_string()));
        all
    }

    fn loaded(file: &str, params: &[&str]) -> OpenTargetsPlugin {
        let mut plugin = OpenTargetsPlugin::new();
        plugin
            .init(&init_params(file, params))
            .expect("the plugin loads");
        plugin
    }

    fn init_error(file: &str, params: &[&str]) -> String {
        let mut plugin = OpenTargetsPlugin::new();
        match plugin.init(&init_params(file, params)) {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected an init error, got {other:?}"),
        }
    }

    fn snv(chr: &str, pos: u64, reference: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.into(),
            pos,
            pos,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn annotated(plugin: &OpenTargetsPlugin, variant: InputVariant) -> InputVariant {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).unwrap();
        variants.remove(0)
    }

    fn keys(headers: &[(String, String)]) -> Vec<&str> {
        headers.iter().map(|(k, _)| k.as_str()).collect()
    }

    /// Without `cols`, the five columns Ensembl VEP's plugin defaults to are
    /// written, each described as a column of the file, the locus-to-gene
    /// score with its meaning first; the variant's lead GWAS record and its
    /// eQTL record print `,`-joined in file order, `NA` where a record has no
    /// value, and the GWAS record that is not a lead variant is left out.
    #[test]
    fn default_columns_are_perl_s() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        let headers = plugin.header_info();
        assert_eq!(
            keys(&headers),
            [
                "OpenTargets_gwasLocusToGeneScore",
                "OpenTargets_gwasGeneId",
                "OpenTargets_gwasDiseases",
                "OpenTargets_qtlGeneId",
                "OpenTargets_qtlBiosampleName",
            ]
        );
        assert_eq!(
            headers[0].1,
            "Locus-to-gene (L2G) scores to predict causal genes at GWAS loci; column from \
             open_targets_vep.tsv.gz"
        );
        assert_eq!(headers[1].1, "column from open_targets_vep.tsv.gz");

        let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
        let data: Vec<(&str, &str)> = variant
            .plugin_data
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        assert_eq!(
            data,
            [
                ("OpenTargets_gwasLocusToGeneScore", "0.5500299334526062,NA"),
                ("OpenTargets_gwasGeneId", "ENSG00000000001,NA"),
                ("OpenTargets_gwasDiseases", "EFO_0004527,NA"),
                ("OpenTargets_qtlGeneId", "NA,ENSG00000000003"),
                ("OpenTargets_qtlBiosampleName", "NA,cerebellum"),
            ]
        );

        let other_allele = annotated(&plugin, snv("21", 5030000, "A", "T"));
        assert!(other_allele.plugin_data.is_empty());
        assert!(other_allele.plugin_records.is_empty());
    }

    /// `cols=all` writes every header column after the four variant columns,
    /// `variantId` to `qtlBiosampleName`, under the `OpenTargets_` prefix.
    #[test]
    fn cols_all_prints_every_header_column() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &["cols=all"]);
        let expected: Vec<String> = HEADER
            .split('\t')
            .skip(4)
            .map(|column| format!("OpenTargets_{column}"))
            .collect();
        assert_eq!(expected.len(), 17);
        assert_eq!(keys(&plugin.header_info()), expected);

        let variant = annotated(&plugin, snv("21", 5030001, "C", "T"));
        assert_eq!(
            variant.plugin_data["OpenTargets_variantId"],
            "21_5030001_C_T"
        );
        assert_eq!(variant.plugin_data["OpenTargets_studyType"], "pqtl");
        assert_eq!(variant.plugin_data["OpenTargets_isLead"], "true");
        assert_eq!(variant.plugin_data["OpenTargets_gwasGeneId"], "NA");
        assert_eq!(variant.plugin_data["OpenTargets_qtlBiosampleName"], "liver");
    }

    /// `OpenTargets.pm:308` and `:318` key the header and the row's values by a
    /// hash over `cols`, so a column named twice is one key: one header line in
    /// its first position and one value.
    #[test]
    fn a_column_named_twice_in_cols_is_declared_once() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(
            &export(dir.path()),
            &["cols=gwasGeneId:qtlBiosampleName:gwasGeneId"],
        );
        assert_eq!(
            keys(&plugin.header_info()),
            ["OpenTargets_gwasGeneId", "OpenTargets_qtlBiosampleName"]
        );
        let variant = annotated(&plugin, snv("21", 5030001, "C", "T"));
        assert_eq!(variant.plugin_data.len(), 2);
        assert_eq!(variant.plugin_data["OpenTargets_qtlBiosampleName"], "liver");
    }

    /// `lead_variants_only=0` keeps the record whose `isLead` is `false`; the
    /// default and `=1` drop it; another value is refused.
    #[test]
    fn lead_variants_only_filters() {
        let dir = tempfile::tempdir().unwrap();
        let file = export(dir.path());
        let all = loaded(&file, &["lead_variants_only=0"]);
        let variant = annotated(&all, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000001,ENSG00000000002,NA"
        );
        assert_eq!(variant.plugin_records["OpenTargets"].len(), 3);

        let leads = loaded(&file, &["lead_variants_only=1"]);
        let variant = annotated(&leads, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000001,NA"
        );

        assert_eq!(
            init_error(&file, &["lead_variants_only=yes"]),
            "OpenTargets: lead_variants_only must be 0 or 1, not 'yes'"
        );
    }

    /// `study_types=0` and `cols=0` are values Perl reads as false
    /// (`if($studytype_request_str)`, OpenTargets.pm:147; `if(!$cols)`,
    /// `:271`), so each is the parameter absent: every family is selected and
    /// the default columns are written, as under `=`.
    #[test]
    fn a_false_study_types_or_cols_value_is_the_parameter_absent() {
        let dir = tempfile::tempdir().unwrap();
        let file = export(dir.path());
        let base = loaded(&file, &[]);
        let base_variant = annotated(&base, snv("21", 5030000, "A", "G"));
        for params in [
            &["study_types=0"][..],
            &["cols=0"][..],
            &["study_types=0", "cols=0"][..],
            &["study_types=", "cols="][..],
        ] {
            let plugin = loaded(&file, params);
            assert_eq!(plugin.header_info(), base.header_info(), "{params:?}");
            let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
            assert_eq!(variant.plugin_data, base_variant.plugin_data, "{params:?}");
        }
    }

    /// `study_types=gwas` keeps the GWAS records and defaults to the GWAS
    /// columns alone, a trailing `:` ignored; `QTL` keeps every QTL type with
    /// the QTL columns; a named QTL type keeps that type; `all` filters
    /// nothing; an unknown name, a named QTL type beside `QTL`, and a family
    /// the file has no columns for are refused.
    #[test]
    fn study_types_filters() {
        let dir = tempfile::tempdir().unwrap();
        let file = export(dir.path());

        let gwas = loaded(&file, &["study_types=gwas"]);
        assert_eq!(
            keys(&gwas.header_info()),
            [
                "OpenTargets_gwasLocusToGeneScore",
                "OpenTargets_gwasGeneId",
                "OpenTargets_gwasDiseases",
            ]
        );
        let variant = annotated(&gwas, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000001"
        );
        assert_eq!(variant.plugin_records["OpenTargets"].len(), 1);
        let trailing_colon = loaded(&file, &["study_types=gwas:"]);
        assert_eq!(keys(&trailing_colon.header_info()).len(), 3);

        let qtl = loaded(&file, &["study_types=QTL"]);
        assert_eq!(
            keys(&qtl.header_info()),
            ["OpenTargets_qtlGeneId", "OpenTargets_qtlBiosampleName"]
        );
        let variant = annotated(&qtl, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlGeneId"],
            "ENSG00000000003"
        );
        let variant = annotated(&qtl, snv("21", 5030001, "C", "T"));
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlGeneId"],
            "ENSG00000000004"
        );

        let eqtl = loaded(&file, &["study_types=eqtl"]);
        let variant = annotated(&eqtl, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlGeneId"],
            "ENSG00000000003"
        );
        let variant = annotated(&eqtl, snv("21", 5030001, "C", "T"));
        assert!(variant.plugin_data.is_empty());

        let all = loaded(&file, &["study_types=all"]);
        assert_eq!(keys(&all.header_info()).len(), 5);
        let variant = annotated(&all, snv("21", 5030000, "A", "G"));
        assert_eq!(variant.plugin_records["OpenTargets"].len(), 2);

        let both = loaded(&file, &["study_types=gwas:eqtl:pqtl"]);
        let variant = annotated(&both, snv("21", 5030001, "C", "T"));
        assert_eq!(variant.plugin_data["OpenTargets_qtlBiosampleName"], "liver");

        assert_eq!(
            init_error(&file, &["study_types=gwas:xyz"]),
            "OpenTargets: unknown study type 'xyz' in 'gwas:xyz'; a study type is all, QTL, \
             gwas or a name ending in qtl"
        );
        assert_eq!(
            init_error(&file, &["study_types=eqtl:QTL"]),
            "OpenTargets: study type 'QTL' in 'eqtl:QTL' conflicts with the QTL study types \
             listed before it"
        );
        assert_eq!(
            init_error(&file, &["study_types=QTL:eqtl"]),
            "OpenTargets: study type 'eqtl' in 'QTL:eqtl' conflicts with the QTL before it, \
             which already selects every QTL study type"
        );

        let gwas_only = format!("#{}", HEADER.rsplit_once("\tqtlGeneId").unwrap().0);
        let mut lines = vec![gwas_only.as_str()];
        lines.extend(RECORDS[..2].iter().map(|r| r.trim_end_matches('\t')));
        let gwas_file = write_tabix_fixture(
            dir.path(),
            "gwas_only.tsv.gz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        assert_eq!(
            init_error(&gwas_file, &["study_types=QTL"]),
            format!(
                "OpenTargets: study type 'QTL' in 'QTL' matches QTL study types but the \
                 header of {gwas_file} has no QTL columns"
            )
        );
        let defaults = loaded(&gwas_file, &[]);
        assert_eq!(keys(&defaults.header_info()).len(), 3);
    }

    /// `QTL` named twice in `study_types` is the one selection it names, every
    /// study type ending in `qtl`, wherever the two stand in the list, so the
    /// plugin loads with the QTL columns and annotates as under `QTL`. Ensembl
    /// VEP's plugin dies on the second token (OpenTargets.pm:161 dereferences
    /// the pattern :170 stored as a list, `Not an ARRAY reference`) and the
    /// plugin is dropped. A named type beside `QTL` stays the conflict both
    /// engines refuse.
    #[test]
    fn a_repeated_qtl_token_selects_the_qtl_study_types_once() {
        let dir = tempfile::tempdir().unwrap();
        let file = export(dir.path());
        for (repeated, once) in [
            ("QTL:QTL", "QTL"),
            ("QTL:gwas:QTL", "QTL:gwas"),
            ("QTL:QTL:QTL", "QTL"),
        ] {
            let base = loaded(&file, &[&format!("study_types={once}")]);
            let plugin = loaded(&file, &[&format!("study_types={repeated}")]);
            assert_eq!(plugin.header_info(), base.header_info(), "{repeated}");
            for site in [snv("21", 5030000, "A", "G"), snv("21", 5030001, "C", "T")] {
                let expected = annotated(&base, site.clone());
                let variant = annotated(&plugin, site);
                assert_eq!(variant.plugin_data, expected.plugin_data, "{repeated}");
                assert_eq!(
                    variant.plugin_records, expected.plugin_records,
                    "{repeated}"
                );
            }
        }
        let qtl = loaded(&file, &["study_types=QTL:QTL"]);
        let variant = annotated(&qtl, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlGeneId"],
            "ENSG00000000003"
        );
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlBiosampleName"],
            "cerebellum"
        );
        assert_eq!(
            init_error(&file, &["study_types=QTL:QTL:eqtl"]),
            "OpenTargets: study type 'eqtl' in 'QTL:QTL:eqtl' conflicts with the QTL before \
             it, which already selects every QTL study type"
        );
    }

    /// The JSON format's records: one object per matched record in file order
    /// under `OpenTargets`, carrying the selected keys as written in the file;
    /// a cell past the line's last filled one is absent (`None`), an empty
    /// cell before a filled one is the empty string.
    #[test]
    fn json_groups_records() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
        let records = &variant.plugin_records["OpenTargets"];
        assert_eq!(records.len(), 2);
        let gwas: Vec<(&str, Option<&str>)> = records[0]
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_deref()))
            .collect();
        assert_eq!(
            gwas,
            [
                (
                    "OpenTargets_gwasLocusToGeneScore",
                    Some("0.5500299334526062")
                ),
                ("OpenTargets_gwasGeneId", Some("ENSG00000000001")),
                ("OpenTargets_gwasDiseases", Some("EFO_0004527")),
                ("OpenTargets_qtlGeneId", None),
                ("OpenTargets_qtlBiosampleName", None),
            ]
        );
        let eqtl: Vec<(&str, Option<&str>)> = records[1]
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_deref()))
            .collect();
        assert_eq!(
            eqtl,
            [
                ("OpenTargets_gwasLocusToGeneScore", Some("")),
                ("OpenTargets_gwasGeneId", Some("")),
                ("OpenTargets_gwasDiseases", Some("")),
                ("OpenTargets_qtlGeneId", Some("ENSG00000000003")),
                ("OpenTargets_qtlBiosampleName", Some("cerebellum")),
            ]
        );
    }

    /// Ensembl VEP splits the line htslib hands it, the line end alone removed,
    /// on tabs (`parse_data`, OpenTargets.pm:396), so the whitespace before the
    /// line end stays on the last cell, and a whitespace-only last cell is a
    /// value: `_join_results` prints it where an empty cell prints `NA`
    /// (`:339`), and the JSON record carries it where an empty cell is `null`.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let header = format!("#{HEADER}");
        let lines = [
            header.as_str(),
            "21\t5030000\tA\tG\t21_5030000_A_G\tgtex_ge_cerebellum\teqtl\tloc3\t6.103\t-7\t1.05262\
             \ttrue\ttrue\t0.95\tSuSie\ttrue\t\t\t\tENSG00000000003\tcerebellum ",
            "21\t5030001\tC\tT\t21_5030001_C_T\tukb_ppp_liver\tpqtl\tloc4\t2.5\t-9\t0.3\ttrue\ttrue\
             \t0.9\tSuSie\ttrue\t\t\t\tENSG00000000004\t ",
        ];
        let file = write_tabix_fixture(
            dir.path(),
            "trailing_space.tsv.gz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let plugin = loaded(&file, &[]);

        let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlBiosampleName"],
            "cerebellum "
        );
        assert_eq!(
            variant.plugin_records["OpenTargets"][0]["OpenTargets_qtlBiosampleName"].as_deref(),
            Some("cerebellum ")
        );

        let variant = annotated(&plugin, snv("21", 5030001, "C", "T"));
        assert_eq!(variant.plugin_data["OpenTargets_qtlBiosampleName"], " ");
        assert_eq!(
            variant.plugin_records["OpenTargets"][0]["OpenTargets_qtlBiosampleName"].as_deref(),
            Some(" ")
        );
    }

    /// A VCF-style deletion record (`AT`/`A` at the anchor base) matches the
    /// variant written as Ensembl VEP stores it (`T`/`-` one base later), and
    /// a VCF-style insertion record its insertion; the records lie one base
    /// before the variant's start and are read by the query window.
    #[test]
    fn vcf_style_indel_records_match_the_trimmed_variant() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        let deletion =
            InputVariant::new("21".into(), 5030011, 5030011, b"T".to_vec(), b"-".to_vec());
        let variant = annotated(&plugin, deletion);
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000005"
        );
        let insertion =
            InputVariant::new("21".into(), 5030021, 5030020, b"-".to_vec(), b"CA".to_vec());
        let variant = annotated(&plugin, insertion);
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000006"
        );
        let other_insertion =
            InputVariant::new("21".into(), 5030021, 5030020, b"-".to_vec(), b"CG".to_vec());
        assert!(annotated(&plugin, other_insertion).plugin_data.is_empty());
    }

    /// Ensembl VEP puts every record its filters keep to the allele match
    /// (`run`, OpenTargets.pm:375-390), which dies on a record whose reference
    /// or alternate cell lies past the line's last filled cell, the undefined
    /// value `split /\t/` leaves (`parse_data`, `:396`; `trim_sequences`,
    /// Utils/Sequence.pm:968). Under `lead_variants_only=0` such a record in
    /// the window two bases before a variant's start to its end fails the
    /// variant wherever it lies among the records that match: `run_batch`
    /// reports the row, the variant carries none of the plugin's fields, and a
    /// variant of the same batch whose window holds no such record keeps its
    /// fields.
    #[test]
    fn a_kept_record_with_no_allele_cell_fails_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let shapes = [
            ("before.tsv.gz", 0, "21\t5030000\tA"),
            ("after.tsv.gz", 3, "21\t5030000\tA"),
            ("trailing_tab.tsv.gz", 0, "21\t5030000\tA\t"),
            ("two_cells.tsv.gz", 0, "21\t5030000"),
        ];
        for (name, index, row) in shapes {
            let file = export_with(dir.path(), name, index, &[row]);
            let plugin = loaded(&file, &["lead_variants_only=0"]);
            let mut variants = vec![
                snv("21", 5030000, "A", "G"),
                snv("21", 5030001, "C", "T"),
                InputVariant::new("21".into(), 5030011, 5030011, b"T".to_vec(), b"-".to_vec()),
            ];
            match plugin.run_batch(&mut variants) {
                Err(PluginError::Run(message)) => assert_eq!(
                    message,
                    format!(
                        "OpenTargets: the row of {file} at 21:5030000 has no reference or \
                         alternate allele cell"
                    ),
                    "{name}"
                ),
                other => panic!("{name}: expected a run error, got {other:?}"),
            }
            for failed in &variants[..2] {
                assert!(failed.plugin_data.is_empty(), "{name}: {}", failed.start);
                assert!(failed.plugin_records.is_empty(), "{name}: {}", failed.start);
            }
            assert_eq!(
                variants[2].plugin_data["OpenTargets_gwasGeneId"], "ENSG00000000005",
                "{name}"
            );
        }
    }

    /// A record the filters drop never reaches the match: under the default
    /// `lead_variants_only=1` a three-cell record has no `isLead` cell and is
    /// dropped (`run`, OpenTargets.pm:352-356), and under a `study_types`
    /// filter it has no `studyType` cell and is dropped (`:359-372`). An empty
    /// cell before a filled one is the defined empty string, which the match
    /// compares as `-` (`trim_sequences`, Utils/Sequence.pm:1033-1034): a lead
    /// record whose alternate cell is empty before its filled cells is no
    /// match for the substitution and is the match for the deletion of its
    /// reference base. In each case the variant keeps the fields of the
    /// records that match and `run_batch` succeeds.
    #[test]
    fn a_filtered_or_empty_allele_record_does_not_fail_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let short = export_with(dir.path(), "short.tsv.gz", 0, &["21\t5030000\tA"]);
        for (params, gwas_gene_ids) in [
            (&[][..], "ENSG00000000001,NA"),
            (&["lead_variants_only=1"][..], "ENSG00000000001,NA"),
            (&["study_types=gwas"][..], "ENSG00000000001"),
            (
                &["lead_variants_only=0", "study_types=gwas"][..],
                "ENSG00000000001,ENSG00000000002",
            ),
        ] {
            let plugin = loaded(&short, params);
            let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
            assert_eq!(
                variant.plugin_data["OpenTargets_gwasGeneId"], gwas_gene_ids,
                "{params:?}"
            );
        }

        let empty_alt = export_with(
            dir.path(),
            "empty_alt.tsv.gz",
            0,
            &[
                "21\t5030000\tA\t\t21_5030000_A_\tGCST9\tgwas\tloc9\t1.0\t-9\t0.1\ttrue\ttrue\
               \t1.0\tPICS\ttrue\tENSG00000000009\t0.9\tEFO_0000009\t\t",
            ],
        );
        for (params, gwas_gene_ids) in [
            (&[][..], "ENSG00000000001,NA"),
            (
                &["lead_variants_only=0"][..],
                "ENSG00000000001,ENSG00000000002,NA",
            ),
        ] {
            let plugin = loaded(&empty_alt, params);
            let substitution = annotated(&plugin, snv("21", 5030000, "A", "G"));
            assert_eq!(
                substitution.plugin_data["OpenTargets_gwasGeneId"], gwas_gene_ids,
                "{params:?}"
            );
            let deletion = annotated(
                &plugin,
                InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"-".to_vec()),
            );
            assert_eq!(
                deletion.plugin_data["OpenTargets_gwasGeneId"], "ENSG00000000009",
                "{params:?}"
            );
            assert_eq!(
                deletion.plugin_records["OpenTargets"].len(),
                1,
                "{params:?}"
            );
        }
    }

    /// A record is matched against every alternate allele of the input record,
    /// so the `T` allele of a record `A/G/T` carries the records its `G` allele
    /// matches, as Ensembl VEP's plugin writes them on every allele's rows.
    #[test]
    fn every_alternate_allele_of_the_record_is_matched() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        let mut other_allele = snv("21", 5030000, "A", "T");
        assert!(annotated(&plugin, other_allele.clone())
            .plugin_data
            .is_empty());
        other_allele.record_allele_string_multi = Some("A/G/T".into());
        let variant = annotated(&plugin, other_allele);
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000001,NA"
        );
    }

    /// A record given on the reverse strand is matched through its
    /// forward-strand alleles, an insertion's included: its inserted bases are
    /// complemented although its `-` reads the same on both strands.
    #[test]
    fn reverse_strand_variant_matches_through_forward_alleles() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        let mut variant = snv("21", 5030001, "G", "A");
        variant.strand = Strand::Reverse;
        let variant = annotated(&plugin, variant);
        assert_eq!(
            variant.plugin_data["OpenTargets_qtlGeneId"],
            "ENSG00000000004"
        );

        let mut insertion =
            InputVariant::new("21".into(), 5030021, 5030020, b"-".to_vec(), b"TG".to_vec());
        insertion.strand = Strand::Reverse;
        let variant = annotated(&plugin, insertion);
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000006"
        );
        let mut as_written =
            InputVariant::new("21".into(), 5030021, 5030020, b"-".to_vec(), b"CA".to_vec());
        as_written.strand = Strand::Reverse;
        assert!(annotated(&plugin, as_written).plugin_data.is_empty());
    }

    /// The export as published carries its header without `#` and an index
    /// that skips the line; the plugin reads that header too.
    #[test]
    fn a_bare_header_file_loads() {
        let dir = tempfile::tempdir().unwrap();
        let mut lines = vec![HEADER];
        lines.extend(RECORDS);
        let file = write_tabix_fixture(
            dir.path(),
            "open_targets_vep_full.tsv.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        );
        let plugin = loaded(&file.display().to_string(), &[]);
        assert_eq!(keys(&plugin.header_info()).len(), 5);
        assert_eq!(
            plugin.header_info()[1].1,
            "column from open_targets_vep_full.tsv.gz"
        );
        let variant = annotated(&plugin, snv("21", 5030000, "A", "G"));
        assert_eq!(
            variant.plugin_data["OpenTargets_gwasGeneId"],
            "ENSG00000000001,NA"
        );
    }

    /// The refusals: no `file`, a positional parameter, a column `cols` names
    /// that the header lacks, a parameter the plugin does not declare.
    #[test]
    fn parameter_refusals_name_the_fault() {
        let dir = tempfile::tempdir().unwrap();
        let file = export(dir.path());
        let mut plugin = OpenTargetsPlugin::new();
        let message = |result: Result<(), PluginError>| match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected an init error, got {other:?}"),
        };
        assert_eq!(
            message(plugin.init(&[])),
            "OpenTargets: no file specified; give file=<path to the Open Targets export>"
        );
        assert_eq!(
            message(plugin.init(std::slice::from_ref(&file))),
            "OpenTargets: no file specified; give file=<path to the Open Targets export>"
        );
        assert_eq!(
            init_error(&file, &["cols=gwasGeneId:score:beta:"]),
            format!("OpenTargets: columns not found in the header of {file}: score")
        );
        assert_eq!(
            init_error(&file, &["cutoff=1"]),
            "OpenTargets: unknown parameter 'cutoff'; OpenTargets accepts file, cols, \
             study_types, lead_variants_only"
        );
        assert_eq!(
            init_error("/no/such/export.tsv.gz", &[]),
            "OpenTargets: data file /no/such/export.tsv.gz not found"
        );
    }

    /// The refusals the header's shape decides: no header line, no column after
    /// the four variant columns, neither a gwas nor a qtl column, no `isLead`
    /// column under `lead_variants_only=1`, no `studyType` column under a
    /// `study_types` filter, and `gwas` against a file with QTL columns alone.
    #[test]
    fn header_refusals_name_the_fault() {
        let dir = tempfile::tempdir().unwrap();
        let shaped = |name: &str, header: Option<&str>| -> String {
            let header = header.map(|h| format!("#{h}"));
            let mut lines: Vec<&str> = header.iter().map(String::as_str).collect();
            lines.extend(RECORDS);
            write_tabix_fixture(dir.path(), name, &lines, IndexSpec::standard())
                .display()
                .to_string()
        };

        let headerless = shaped("headerless.tsv.gz", None);
        assert_eq!(
            init_error(&headerless, &[]),
            format!("OpenTargets: no header line found in {headerless}")
        );

        let variant_columns = HEADER.split('\t').take(4).collect::<Vec<_>>().join("\t");
        let bare = shaped("variant_columns.tsv.gz", Some(&variant_columns));
        assert_eq!(
            init_error(&bare, &[]),
            format!(
                "OpenTargets: the header of {bare} has no columns after the four variant columns"
            )
        );

        let no_family = shaped(
            "no_family.tsv.gz",
            Some(&format!("{variant_columns}\tvariantId")),
        );
        assert_eq!(
            init_error(&no_family, &[]),
            format!(
                "OpenTargets: {no_family} carries neither GWAS nor QTL data: its header has no \
                 gwas or qtl columns"
            )
        );

        let no_is_lead = shaped("no_is_lead.tsv.gz", Some(&HEADER.replace("\tisLead", "")));
        assert_eq!(
            init_error(&no_is_lead, &[]),
            format!(
                "OpenTargets: lead_variants_only=1 needs an isLead column, which the header of \
                 {no_is_lead} lacks; pass lead_variants_only=0"
            )
        );
        let unfiltered = loaded(&no_is_lead, &["lead_variants_only=0"]);
        assert_eq!(keys(&unfiltered.header_info()).len(), 5);

        let no_study_type = shaped(
            "no_study_type.tsv.gz",
            Some(&HEADER.replace("\tstudyType", "")),
        );
        assert_eq!(
            init_error(&no_study_type, &["study_types=gwas"]),
            format!(
                "OpenTargets: study_types needs a studyType column, which the header of \
                 {no_study_type} lacks"
            )
        );

        let qtl_only = shaped(
            "qtl_only.tsv.gz",
            Some(&HEADER.replace("\tgwasGeneId\tgwasLocusToGeneScore\tgwasDiseases", "")),
        );
        assert_eq!(
            init_error(&qtl_only, &["study_types=gwas"]),
            format!(
                "OpenTargets: study type 'gwas' in 'gwas' matches GWAS study types but the \
                 header of {qtl_only} has no GWAS columns"
            )
        );
    }

    /// The trimmed forms of a pair of alleles, as `trim_sequences` leaves them:
    /// the shared bases gone from both ends, the start past the left ones, `-`
    /// for an allele with nothing left; a one-base pair has one form.
    #[test]
    fn trimming_follows_trim_sequences() {
        let at = KeyPosition::Number;
        assert_eq!(
            trimmed_forms(b"A", b"G", at(100)),
            [(b"A".to_vec(), b"G".to_vec(), at(100))]
        );
        assert_eq!(
            trimmed_forms(b"ATG", b"AAG", at(100)),
            [
                (b"T".to_vec(), b"A".to_vec(), at(101)),
                (b"T".to_vec(), b"A".to_vec(), at(101)),
            ]
        );
        assert_eq!(
            trimmed_forms(b"AT", b"A", at(100)),
            [
                (b"T".to_vec(), b"-".to_vec(), at(101)),
                (b"T".to_vec(), b"-".to_vec(), at(101)),
            ]
        );
        // `TA`/`TATA`: left first leaves `-`/`TA` at 102, right first `-`/`TA` at 100.
        assert_eq!(
            trimmed_forms(b"TA", b"TATA", at(100)),
            [
                (b"-".to_vec(), b"TA".to_vec(), at(102)),
                (b"-".to_vec(), b"TA".to_vec(), at(100)),
            ]
        );
    }

    /// Perl's `split` drops the empty pieces at the end of a list and keeps
    /// the ones before a filled piece; a list of separators alone, or an empty
    /// list, has no pieces.
    #[test]
    fn perl_split_drops_trailing_empty_pieces_only() {
        assert_eq!(perl_split("a::b", ':'), ["a", "", "b"]);
        assert_eq!(perl_split(":a", ':'), ["", "a"]);
        assert_eq!(perl_split("a:", ':'), ["a"]);
        assert!(perl_split(":", ':').is_empty());
        assert!(perl_split("", ':').is_empty());
    }

    /// Ensembl VEP 116.2 never hands OpenTargets a structural variant: the
    /// module keeps the default `variant_feature_types` of
    /// `['VariationFeature']` (BaseVepPlugin.pm:100), which the dispatch checks
    /// before running a plugin on a line (OutputFactory.pm:2220). A `<DEL>`
    /// over the export's records gets no field, and nor does a matching
    /// substitution carrying the structural mark, which pins the gate rather
    /// than the allele match; the plain substitution is annotated.
    #[test]
    fn a_structural_variant_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = loaded(&export(dir.path()), &[]);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = snv("21", 5030000, "N", "<DEL>");
        deletion.end = 5030100;
        deletion.is_structural = true;
        let deletion = annotated(&plugin, deletion);
        assert!(deletion.plugin_data.is_empty());
        assert!(deletion.plugin_records.is_empty());
        let mut flagged = snv("21", 5030000, "A", "G");
        flagged.is_structural = true;
        assert!(annotated(&plugin, flagged).plugin_data.is_empty());
        assert!(!annotated(&plugin, snv("21", 5030000, "A", "G"))
            .plugin_data
            .is_empty());
    }

    /// A record's position cell is keyed by `get_matched_variant_alleles`
    /// (`OpenTargets.pm:376-388`) as text when the alleles lose no leading base
    /// (ensembl-variation Utils/Sequence.pm:1232), and as the text `$start++`
    /// (`:405`) steps a cell no arithmetic has read to once `parse_data` drops
    /// a shared leading base: for a substitution `6000000` annotates and
    /// `6000000.0` and `6000000.5` do not, nor does `6000001x` at 6000001; for
    /// the deletion `AT>A` at `6000010.0` and `6000030x` the stepped start is
    /// the variant's, while `6000020.5` steps to 6000021.5 and `06000040` to
    /// the text `06000041`, and each matches nothing.
    #[test]
    fn a_position_cell_is_matched_as_the_key_holds_it() {
        let dir = tempfile::tempdir().unwrap();
        let gwas = |pos: &str, reference: &str, alt: &str, gene: &str| {
            format!(
                "21\t{pos}\t{reference}\t{alt}\tid\tGCST9\tgwas\tloc9\t3.0\t-8\t0.1\ttrue\ttrue\t1.0\
                 \tPICS\ttrue\t{gene}\t0.7\tEFO_0000009\t\t"
            )
        };
        let rows = [
            gwas("6000000", "A", "G", "ENSG00000000010"),
            gwas("6000000.0", "A", "C", "ENSG00000000011"),
            gwas("6000000.5", "A", "T", "ENSG00000000012"),
            gwas("6000001x", "A", "G", "ENSG00000000013"),
            gwas("6000010.0", "AT", "A", "ENSG00000000014"),
            gwas("6000020.5", "AT", "A", "ENSG00000000015"),
            gwas("6000030x", "AT", "A", "ENSG00000000016"),
            gwas("06000040", "AT", "A", "ENSG00000000017"),
        ];
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let file = export_with(dir.path(), "cells.tsv.gz", RECORDS.len(), &rows);
        let plugin = loaded(&file, &[]);
        let gene = |variant: InputVariant| {
            annotated(&plugin, variant)
                .plugin_data
                .get("OpenTargets_gwasGeneId")
                .cloned()
        };
        assert_eq!(
            gene(snv("21", 6000000, "A", "G")).as_deref(),
            Some("ENSG00000000010")
        );
        assert_eq!(gene(snv("21", 6000000, "A", "C")), None);
        assert_eq!(gene(snv("21", 6000000, "A", "T")), None);
        assert_eq!(gene(snv("21", 6000001, "A", "G")), None);
        let deletion =
            |start: u64| InputVariant::new("21".into(), start, start, b"T".to_vec(), b"-".to_vec());
        assert_eq!(gene(deletion(6000011)).as_deref(), Some("ENSG00000000014"));
        assert_eq!(gene(deletion(6000021)), None);
        assert_eq!(gene(deletion(6000031)).as_deref(), Some("ENSG00000000016"));
        assert_eq!(
            gene(deletion(6000041)),
            None,
            "06000040 steps to the text 06000041"
        );

        // A Perl-false cell returns the match before an undefined allele would
        // fail it (Sequence.pm:1173 before :1230); a Perl-true one still fails.
        let record = |pos: &str, cells: &[&str]| TabixRecord {
            chr: "21".to_string(),
            start: 1,
            end: 1,
            columns: ["21", pos]
                .into_iter()
                .chain(cells.iter().copied())
                .map(str::to_string)
                .collect(),
            header: None,
        };
        let variant = snv("21", 1, "A", "G");
        let alternates = record_alternate_alleles(&variant);
        assert_eq!(
            alleles_match(&variant, &alternates, &record("0", &["A"])),
            Some(false)
        );
        assert_eq!(
            alleles_match(&variant, &alternates, &record("", &["A"])),
            Some(false)
        );
        assert_eq!(
            alleles_match(&variant, &alternates, &record("1", &["A"])),
            None
        );
        let deletion = InputVariant::new("21".into(), 1, 1, b"A".to_vec(), b"-".to_vec());
        let alternates = record_alternate_alleles(&deletion);
        assert_eq!(
            alleles_match(&deletion, &alternates, &record("0", &["GA", "G"])),
            Some(true)
        );
        // `abc` steps to the text `abd`, where a numified cell would reach 1;
        // `-1` steps to 0, Perl-false, and ends the match.
        assert_eq!(
            alleles_match(&deletion, &alternates, &record("abc", &["GA", "G"])),
            Some(false)
        );
        assert_eq!(
            alleles_match(&deletion, &alternates, &record("-1", &["GA", "G"])),
            Some(false)
        );
    }
}
