// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! gnomADMt plugin: gnomAD mitochondrial allele frequencies from the gnomAD
//! chrM sites VCF, on every row of a matching record.
//!
//! The output fields are INFO fields of the data file, named `gnomAD_<field>`
//! and described by the file's own `##INFO` descriptions (gnomADMt.pm:146-163,
//! :217-249); `fields=` selects them, `all` every one. `hap_filter=` and
//! `pop_filter=` keep the haplogroups or populations named, in the order
//! named, or every one in file order under `all` (`:82`), of every
//! per-haplogroup or per-population field, whose values are `|`-joined lists
//! in the order the `hap_AN` or `pop_AN` description states (`:104-130`,
//! `:165-205`, `:340-364`). The plugin declares `Feature` and
//! `Intergenic` rows (`:214`), so a record's values print on every row of it.

use std::io::BufReader;
use std::path::Path;

use indexmap::IndexMap;
use noodles::bgzf;
use tracing::{debug, warn};

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, lossy_text, BatchQueryResult, BoundedLines, TabixAnnotator, TabixConfig,
    TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;
use vep_core::vcf_header::fileformat_refusal;

const NAME: &str = "gnomADMt";

/// The prefix of every output field (`:207`).
const FIELD_PREFIX: &str = "gnomAD";

/// The fields written without `fields=` (`:159`).
const DEFAULT_FIELDS: &str = "hap_AC_hom:hap_AC_het:hap_AF_hom:hap_AF_het:hap_AN";

/// The VCF columns the plugin reads (`:304`).
const POS_COL: usize = 1;
const REF_COL: usize = 3;
const ALT_COL: usize = 4;
const INFO_COL: usize = 7;

/// Bases read before the record's start (`:267`).
const WINDOW_BEFORE: u64 = 2;

/// A group of per-haplogroup or per-population fields: the prefix its field
/// names open with, the field whose description states the group order, the
/// group noun (`<noun> order:` in that description, singular and plural in
/// the messages), the field whose list is `,`-joined where the others are
/// `|`-joined (`:344`, `:357`), and the `*_filter=` parameter.
struct GroupKind {
    prefix: &'static str,
    order_field: &'static str,
    item: &'static str,
    items: &'static str,
    comma_field: &'static str,
    parameter: &'static str,
}

const HAPLOGROUPS: GroupKind = GroupKind {
    prefix: "hap_",
    order_field: "hap_AN",
    item: "haplogroup",
    items: "haplogroups",
    comma_field: "hap_hl_hist",
    parameter: "hap_filter",
};

const POPULATIONS: GroupKind = GroupKind {
    prefix: "pop_",
    order_field: "pop_AN",
    item: "population",
    items: "populations",
    comma_field: "pop_hl_hist",
    parameter: "pop_filter",
};

/// One selected INFO field: its id, its output key and the description the
/// headers print, the file's with a filtered group's order rewritten.
#[derive(Debug, Clone)]
struct Field {
    id: String,
    key: String,
    description: String,
}

/// The positions, in a field's list, of the groups a filter keeps, in the order
/// the filter names them (`:176-182`, `:197-203`).
#[derive(Debug, Clone, Default)]
struct GroupFilter {
    indexes: Vec<usize>,
}

/// gnomADMt plugin.
pub struct GnomadMtPlugin {
    store: Option<Box<dyn AnnotationStore>>,
    /// The data file, as the failure on one of its rows names it.
    file: String,
    /// The selected fields in the order given, each once.
    fields: Vec<Field>,
    haplogroups: Option<GroupFilter>,
    populations: Option<GroupFilter>,
}

impl Default for GnomadMtPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl GnomadMtPlugin {
    /// Create an uninitialized plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            file: String::new(),
            fields: Vec::new(),
            haplogroups: None,
            populations: None,
        }
    }

    /// The region read for a variant (`:263-267`): two bases before the start
    /// to the end, the end raised to the start for an insertion, whose end
    /// precedes its start (`:265`).
    fn window(variant: &InputVariant) -> (u64, u64) {
        let high = variant.start.max(variant.end);
        (variant.start.saturating_sub(WINDOW_BEFORE).max(1), high)
    }

    /// Resolves a `*_filter=` parameter against the group order the file
    /// states, or `None` when no selected field is of the group, which Ensembl
    /// VEP warns about and ignores (`:168-170`, `:189-191`).
    fn resolve_filter(
        &mut self,
        kind: &GroupKind,
        requested: &str,
        info: &[(String, String)],
        file: &str,
    ) -> Result<Option<GroupFilter>, PluginError> {
        if !self
            .fields
            .iter()
            .any(|f| is_group_field(&f.id, kind.prefix))
        {
            warn!(
                "{NAME}: {} filtering can only be applied to {} fields.",
                capitalized(kind.item),
                kind.item
            );
            return Ok(None);
        }
        let available = info
            .iter()
            .find(|(id, _)| id == kind.order_field)
            .and_then(|(_, description)| group_order(description, kind.item))
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: {} needs the {} order list in the {} INFO description of {file}",
                    kind.parameter, kind.item, kind.order_field
                ))
            })?;
        // `:82`: `all` is every group the file states, in file order, before
        // any name is validated.
        let valid: Vec<&str> = if requested == "all" {
            available.iter().map(String::as_str).collect()
        } else {
            let (valid, invalid): (Vec<&str>, Vec<&str>) = requested
                .split(':')
                .partition(|name| available.iter().any(|a| a == name));
            if valid.is_empty() {
                return Err(PluginError::Init(format!(
                    "{NAME}: all {} given are invalid. Available {} are: {}",
                    kind.items,
                    kind.items,
                    available.join(", ")
                )));
            }
            if !invalid.is_empty() {
                warn!(
                    "{NAME}: the following {} are not valid and were ignored: {}",
                    kind.item,
                    invalid.join(", ")
                );
            }
            valid
        };
        let indexes = valid
            .iter()
            .filter_map(|name| available.iter().position(|a| a == name))
            .collect();
        let rewritten = format!(
            "{} order: [{}]",
            kind.item,
            valid
                .iter()
                .map(|name| format!("'{name}'"))
                .collect::<Vec<_>>()
                .join(",")
        );
        for field in self
            .fields
            .iter_mut()
            .filter(|f| is_group_field(&f.id, kind.prefix))
        {
            field.description = rewrite_group_order(&field.description, kind.item, &rewritten);
        }
        Ok(Some(GroupFilter { indexes }))
    }

    /// The value of `field` on a row's INFO column, its group list filtered and
    /// reordered when a filter applies (`:340-364`); the empty value for a field
    /// the row lacks or carries as a flag, which Ensembl VEP leaves undefined.
    fn field_value(&self, field: &Field, info: &[(&str, &str)]) -> String {
        let Some((_, value)) = info.iter().find(|(key, _)| *key == field.id) else {
            return String::new();
        };
        if value.is_empty() {
            return String::new();
        }
        for (kind, filter) in [
            (&HAPLOGROUPS, &self.haplogroups),
            (&POPULATIONS, &self.populations),
        ] {
            let Some(filter) = filter else {
                continue;
            };
            if !field.id.starts_with(kind.prefix) {
                continue;
            }
            let separator = if field.id == kind.comma_field {
                ","
            } else {
                "|"
            };
            let values: Vec<&str> = value.split(separator).collect();
            return filter
                .indexes
                .iter()
                .map(|&i| values.get(i).copied().unwrap_or(""))
                .collect::<Vec<_>>()
                .join(separator);
        }
        value.to_string()
    }

    /// The selected fields' values of the first row of `records`, read in file
    /// order, whose alleles match the variant's, both sides trimmed of their
    /// shared bases (`:275-296`), through `allele_match::matched_variant_alleles`
    /// on the row as `_parse_vcf` leaves it; `Ok(None)` when none does. A row
    /// read before that match with no reference or alternate cell is the row
    /// Ensembl VEP's match dies on (`trim_sequences`, ensembl-variation
    /// Utils/Sequence.pm:968, reached from `:1230`), so it fails the variant,
    /// whose line then carries none of the fields; the loop stops at its first
    /// match (`last`, `:294`), so such a row behind it is never read. `:304`
    /// assigns `split /\t/` to a list of eight scalars, a split with a limit
    /// that keeps trailing empty cells, so a cell is missing only when the row
    /// ends before it: an empty cell, trailing or not, is the defined empty
    /// string and matches as `-`.
    fn matching_values(
        &self,
        records: &[&TabixRecord],
        variant: &InputVariant,
    ) -> Result<Option<Vec<(String, String)>>, PluginError> {
        let mut matched = None;
        for record in records {
            let undefined_allele =
                record.columns.get(REF_COL).is_none() || record.columns.get(ALT_COL).is_none();
            // A Perl-false position cell, which an undefined allele never moves,
            // returns the match before the allele dies (Sequence.pm:1173 before
            // :1230).
            if undefined_allele && !allele_match::perl_false_position(record, POS_COL) {
                return Err(allele_match::undefined_allele_row(NAME, &self.file, record));
            }
            if allele_match::matched_variant_alleles(
                &allele_match::vcf_adjusted_row(record, POS_COL, REF_COL, ALT_COL),
                variant,
                POS_COL,
                REF_COL,
                ALT_COL,
            ) {
                matched = Some(*record);
                break;
            }
        }
        let Some(record) = matched else {
            return Ok(None);
        };
        let info: Vec<(&str, &str)> = record
            .columns
            .get(INFO_COL)
            .map(|column| column.split(';').map(split_info_pair).collect())
            .unwrap_or_default();
        Ok(Some(
            self.fields
                .iter()
                .map(|field| (field.key.clone(), self.field_value(field, &info)))
                .collect(),
        ))
    }
}

/// Whether `id` is a field of the group opening with `prefix`: it holds the
/// prefix followed by at least one character (`/hap_.+/`, `:168`, `:229`).
fn is_group_field(id: &str, prefix: &str) -> bool {
    id.split_once(prefix)
        .is_some_and(|(_, after_prefix)| !after_prefix.is_empty())
}

/// `word` with its first letter in upper case.
fn capitalized(word: &str) -> String {
    let mut chars = word.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().chain(chars).collect(),
        None => String::new(),
    }
}

/// One `key=value` of an INFO column; a flag is the key with the empty value.
fn split_info_pair(pair: &str) -> (&str, &str) {
    pair.split_once('=').unwrap_or((pair, ""))
}

/// The group names a description states as `<word> order: ['A', 'B']`
/// (`:109-111`): the quoted words inside the brackets, `None` without the
/// statement.
fn group_order(description: &str, word: &str) -> Option<Vec<String>> {
    let statement = format!("{word} order: [");
    let start = description.find(&statement)? + statement.len();
    let end = start + description[start..].find(']')?;
    let names: Vec<String> = description[start..end]
        .split('\'')
        .skip(1)
        .step_by(2)
        .filter(|name| !name.is_empty() && name.chars().all(|c| c.is_alphanumeric() || c == '_'))
        .map(str::to_string)
        .collect();
    Some(names)
}

/// `description` with its `<word> order: [...]` statement replaced by
/// `rewritten` (`:125-127`); unchanged without the statement.
fn rewrite_group_order(description: &str, word: &str, rewritten: &str) -> String {
    let statement = format!("{word} order: [");
    let Some(start) = description.find(&statement) else {
        return description.to_string();
    };
    let Some(close) = description[start..].find(']') else {
        return description.to_string();
    };
    format!(
        "{}{rewritten}{}",
        &description[..start],
        &description[start + close + 1..]
    )
}

/// The `(ID, Description)` of every `##INFO` line of a bgzip-compressed VCF,
/// in file order; a line without a quoted description describes its field as
/// the empty text. The plugin opens its file through ensembl-io's VCF parser
/// (`gnomADMt.pm:146`), which ends the run on a `##fileformat` line naming a
/// major version other than 4 or no `<digits>.<digits>` version
/// (`BaseVCF4.pm:125, 129`), so such a line refuses the file with Ensembl VEP's
/// text behind the line read; each line is read less its `\n` or `\r\n`, as
/// htslib hands it over.
fn info_fields(path: &Path) -> Result<Vec<(String, String)>, PluginError> {
    let file = std::fs::File::open(path).map_err(|e| {
        PluginError::Init(format!("{NAME}: failed to open {}: {e}", path.display()))
    })?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));
    let mut fields = Vec::new();
    for line in BoundedLines::new(reader, path) {
        let line = line?;
        if !line.starts_with(b"#") {
            break;
        }
        if let Some(refusal) = fileformat_refusal(&line) {
            return Err(PluginError::Init(format!(
                "{NAME}: {}: {}: {refusal}",
                path.display(),
                lossy_text(line)
            )));
        }
        let line = lossy_text(line);
        let Some(body) = line
            .strip_prefix("##INFO=<")
            .and_then(|rest| rest.strip_suffix('>'))
        else {
            continue;
        };
        let mut id = None;
        let mut description = String::new();
        for (key, value) in header_pairs(body) {
            match key {
                "ID" => id = Some(value),
                "Description" => description = value,
                _ => {}
            }
        }
        if let Some(id) = id {
            fields.push((id, description));
        }
    }
    Ok(fields)
}

/// The `key=value` pairs of a VCF header line's `<...>` body: a value is
/// quoted, holding any character but an unescaped `"`, or runs to the next
/// comma.
fn header_pairs(body: &str) -> Vec<(&str, String)> {
    let mut pairs = Vec::new();
    let mut rest = body;
    while !rest.is_empty() {
        let Some((key, after)) = rest.split_once('=') else {
            break;
        };
        if let Some(quoted) = after.strip_prefix('"') {
            let mut value = String::new();
            let mut chars = quoted.char_indices();
            let mut consumed = quoted.len();
            while let Some((i, c)) = chars.next() {
                match c {
                    '\\' => {
                        if let Some((_, escaped)) = chars.next() {
                            value.push(escaped);
                        }
                    }
                    '"' => {
                        consumed = i + 1;
                        break;
                    }
                    other => value.push(other),
                }
            }
            pairs.push((key, value));
            rest = quoted[consumed..].trim_start_matches(',');
        } else {
            let (value, remainder) = after.split_once(',').unwrap_or((after, ""));
            pairs.push((key, value.to_string()));
            rest = remainder;
        }
    }
    pairs
}

impl BuiltinPlugin for GnomadMtPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`gnomADMt.pm:213-214`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        &RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.fields
            .iter()
            .map(|field| (field.key.clone(), field.description.clone()))
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(
            NAME,
            &["file", "fields", "hap_filter", "pop_filter"],
            params,
        )?;
        if !params.positional().is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: parameters are key=value: file=path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz[,fields=a:b|all][,hap_filter=HV:A][,pop_filter=eas:fin]"
            )));
        }
        let Some(file) = params.get("file").filter(|f| !f.is_empty()) else {
            return Err(PluginError::Init(format!(
                "{NAME}: no file specified; give file=path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
            )));
        };
        if !Path::new(file).is_file() {
            return Err(PluginError::Init(format!(
                "{NAME}: data file {file} not found"
            )));
        }
        // gnomADMt.pm:146 reads the header through tabix, which needs an index,
        // before `add_file` checks for `<file>.tbi` (:156); `info_fields` needs
        // none, so the check comes first and an unindexed file is refused for
        // its index as there.
        check_table_or_store(Path::new(file))?;

        let info = info_fields(Path::new(file))?;
        let has = |id: &str| info.iter().any(|(field, _)| field == id);
        if !has("AF_hom") || !has("AF_het") {
            return Err(PluginError::Init(format!(
                "{NAME}: {file} does not carry the AF_hom and AF_het INFO fields; give the file as gnomAD distributes it"
            )));
        }

        let requested = params.get("fields").unwrap_or(DEFAULT_FIELDS);
        let selected: Vec<&str> = if requested == "all" {
            info.iter().map(|(id, _)| id.as_str()).collect()
        } else {
            let (valid, invalid): (Vec<&str>, Vec<&str>) =
                requested.split(':').partition(|name| has(name));
            if valid.is_empty() {
                return Err(PluginError::Init(format!(
                    "{NAME}: all fields given are invalid. Available fields are: {}",
                    info.iter()
                        .map(|(id, _)| id.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )));
            }
            if !invalid.is_empty() {
                warn!(
                    "{NAME}: the following fields are not valid and were ignored: {}",
                    invalid.join(", ")
                );
            }
            valid
        };
        self.fields.clear();
        for id in selected {
            if self.fields.iter().any(|f| f.id == id) {
                continue;
            }
            let description = info
                .iter()
                .find(|(field, _)| field == id)
                .map(|(_, description)| description.clone())
                .unwrap_or_default();
            self.fields.push(Field {
                id: id.to_string(),
                key: format!("{FIELD_PREFIX}_{id}"),
                description,
            });
        }

        self.haplogroups = match params.get("hap_filter") {
            Some(requested) => self.resolve_filter(&HAPLOGROUPS, requested, &info, file)?,
            None => None,
        };
        self.populations = match params.get("pop_filter") {
            Some(requested) => self.resolve_filter(&POPULATIONS, requested, &info, file)?,
            None => None,
        };

        let config = TabixConfig {
            file_path: file.into(),
            chr_col: 0,
            start_col: POS_COL,
            end_col: Some(INFO_COL),
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        // A row keeps the whitespace before its line end, which Ensembl VEP's
        // `split /\t/` leaves on the INFO column and its `split /=/` on the
        // value of the column's last pair (`:304`, `:327-329`).
        // Opened on no assembly, so a mitochondrial record reads the table under
        // every spelling: the module itself names such a record `M` to its table
        // (gnomADMt.pm:256-261).
        self.store = Some(open_best_store_with(
            Path::new(file),
            config,
            None,
            TabixAnnotator::keeping_trailing_whitespace,
        )?);
        self.file = file.to_string();
        debug!(
            file,
            fields = self.fields.len(),
            "gnomADMt plugin initialized"
        );
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
            .map(|v| {
                let (low, high) = Self::window(v);
                (v.chr.as_str(), low, high)
            })
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
            let (low, high) = Self::window(variant);
            let records = batch.lookup(&variant.chr, low, high);
            if records.is_empty() {
                continue;
            }
            match self.matching_values(&records, variant) {
                Ok(Some(values)) => variant.plugin_data.extend(values),
                Ok(None) => {}
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
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_bgzf_bytes, write_tabix_fixture, IndexSpec};
    use std::path::PathBuf;
    use vep_core::consequence::FeatureType;

    const HAP_ORDER: &str = "haplogroup order: ['A', 'B', 'HV']";
    const POP_ORDER: &str = "population order: ['afr', 'eas', 'fin']";

    /// A gnomAD chrM sites VCF of the real header shape: three haplogroups,
    /// three populations, a flag field, and rows on `chrM`.
    fn rows() -> Vec<String> {
        let mut lines = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=filters,Number=.,Type=String,Description=\"Site-level filters\">".to_string(),
            "##INFO=<ID=hap_defining_variant,Number=0,Type=Flag,Description=\"Present if variant is present as a haplogroup defining variant\">".to_string(),
            "##INFO=<ID=AC_hom,Number=1,Type=Integer,Description=\"Allele count restricted to variants with a heteroplasmy level >= 0.95\">".to_string(),
            "##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"Allele frequency restricted to variants with a heteroplasmy level >= 0.95\">".to_string(),
            "##INFO=<ID=AF_het,Number=1,Type=Float,Description=\"Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < 0.95\">".to_string(),
            format!("##INFO=<ID=hap_AN,Number=1,Type=String,Description=\"List of overall allele number for each haplogroup, {HAP_ORDER}\">"),
            format!("##INFO=<ID=hap_AC_het,Number=1,Type=String,Description=\"List of AC_het for each haplogroup, {HAP_ORDER}\">"),
            format!("##INFO=<ID=hap_AC_hom,Number=1,Type=String,Description=\"List of AC_hom for each haplogroup, {HAP_ORDER}\">"),
            format!("##INFO=<ID=hap_AF_hom,Number=1,Type=String,Description=\"List of AF_hom for each haplogroup, {HAP_ORDER}\">"),
            format!("##INFO=<ID=hap_AF_het,Number=1,Type=String,Description=\"List of AF_het for each haplogroup, {HAP_ORDER}\">"),
            format!("##INFO=<ID=hap_hl_hist,Number=3,Type=String,Description=\"Histogram of heteroplasmy levels for each haplogroup; bin edges are: [0.0, 0.5, 1.0], {HAP_ORDER}\">"),
            format!("##INFO=<ID=pop_AN,Number=1,Type=String,Description=\"List of overall allele number for each population, {POP_ORDER}\">"),
            format!("##INFO=<ID=pop_AC_hom,Number=1,Type=String,Description=\"List of AC_hom for each population, {POP_ORDER}\">"),
            format!("##INFO=<ID=pop_hl_hist,Number=3,Type=String,Description=\"Histogram of heteroplasmy levels for each population; bin edges are: [0.0, 0.5, 1.0], {POP_ORDER}\">"),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let info = |ac_hom: u32, filters: &str, flag: bool| {
            format!(
                "{}{}AC_hom={ac_hom};AF_hom=0.5;AF_het=0.25;hap_AN=10|20|30;hap_AC_het=1|2|3;hap_AC_hom=4|5|6;hap_AF_hom=0.4|0.5|0.6;hap_AF_het=0.1|0.2|0.3;hap_hl_hist=1|0|0,0|2|0,0|0|3;pop_AN=100|200|300;pop_AC_hom=7|8|9;pop_hl_hist=1|1|1,2|2|2,3|3|3",
                if filters.is_empty() { String::new() } else { format!("filters={filters};") },
                if flag { "hap_defining_variant;" } else { "" }
            )
        };
        lines.push(format!("chrM\t3\t.\tT\tC\t.\tPASS\t{}", info(19, "", true)));
        lines.push(format!(
            "chrM\t21\t.\tAT\tA\t.\tnpg\t{}",
            info(2, "npg", false)
        ));
        lines.push(format!(
            "chrM\t42\trs377245343\tT\tTC\t.\tPASS\t{}",
            info(88, "", false)
        ));
        lines.push(format!(
            "chrM\t500\t.\tAA\tA\t.\tPASS\t{}",
            info(7, "", false)
        ));
        lines.push(format!(
            "chrM\t600\t.\tAA\tAAA\t.\tPASS\t{}",
            info(5, "", false)
        ));
        lines.push(format!(
            "chrM\t8993\t.\tT\tC\t.\tPASS\t{}",
            info(1, "", false)
        ));
        lines.push(format!(
            "chrM\t8993\t.\tT\tG\t.\tnpg\t{}",
            info(0, "npg", false)
        ));
        lines
    }

    fn fixture(dir: &Path) -> PathBuf {
        let lines = rows();
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        write_tabix_fixture(dir, "chrM.vcf.bgz", &refs, IndexSpec::standard())
    }

    /// The fixture's rows with `inserted` merged in by position, each ahead of
    /// the rows at its position when `ahead` is set and behind them otherwise,
    /// written to `file_name`.
    fn fixture_with_rows(dir: &Path, file_name: &str, inserted: &[&str], ahead: bool) -> PathBuf {
        let position = |line: &str| -> u64 { line.split('\t').nth(1).unwrap().parse().unwrap() };
        let mut lines = rows();
        for row in inserted {
            let at = lines
                .iter()
                .position(|line| {
                    !line.starts_with('#')
                        && (position(line) > position(row)
                            || (ahead && position(line) == position(row)))
                })
                .unwrap_or(lines.len());
            lines.insert(at, row.to_string());
        }
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        write_tabix_fixture(dir, file_name, &refs, IndexSpec::standard())
    }

    /// The plugin on `file` selecting the `AC_hom` field alone.
    fn init_file(file: &Path) -> GnomadMtPlugin {
        let mut plugin = GnomadMtPlugin::new();
        plugin
            .init(&[
                format!("file={}", file.display()),
                "fields=AC_hom".to_string(),
            ])
            .unwrap();
        plugin
    }

    /// Runs the batch expecting the failure of a row with no allele cell: the
    /// `Run` error's text names it, and the variants come back as annotated.
    fn annotate_failing_on_an_undefined_row(
        plugin: &GnomadMtPlugin,
        mut variants: Vec<InputVariant>,
    ) -> Vec<InputVariant> {
        match plugin.run_batch(&mut variants) {
            Err(PluginError::Run(message)) => assert!(
                message.contains("has no reference or alternate allele cell"),
                "{message}"
            ),
            other => panic!("expected PluginError::Run, got {other:?}"),
        }
        variants
    }

    fn init_with(dir: &Path, extra: &[&str]) -> Result<GnomadMtPlugin, PluginError> {
        let file = fixture(dir);
        let mut params = vec![format!("file={}", file.display())];
        params.extend(extra.iter().map(|p| p.to_string()));
        let mut plugin = GnomadMtPlugin::new();
        plugin.init(&params)?;
        Ok(plugin)
    }

    /// The `Init` message `result` failed with; any other outcome fails the test.
    fn expect_init_error<T>(result: Result<T, PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            Ok(_) => panic!("expected PluginError::Init"),
            Err(other) => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// The `Init` message of initializing the fixture with `extra` parameters.
    fn init_message(dir: &Path, extra: &[&str]) -> String {
        expect_init_error(init_with(dir, extra))
    }

    fn keys(plugin: &GnomadMtPlugin) -> Vec<String> {
        plugin.header_info().into_iter().map(|(k, _)| k).collect()
    }

    fn description(plugin: &GnomadMtPlugin, key: &str) -> String {
        plugin
            .header_info()
            .into_iter()
            .find(|(k, _)| k == key)
            .map(|(_, d)| d)
            .unwrap_or_else(|| panic!("no header field {key}"))
    }

    fn variant(chr: &str, start: u64, end: u64, reference: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.into(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn annotate(plugin: &GnomadMtPlugin, mut variants: Vec<InputVariant>) -> Vec<InputVariant> {
        plugin.run_batch(&mut variants).unwrap();
        variants
    }

    fn value(variant: &InputVariant, key: &str) -> Option<String> {
        variant.plugin_data.get(key).cloned()
    }

    /// A file lacking both its index and its INFO fields is refused naming the
    /// index, as gnomADMt.pm refuses it: the header read (:146) needs an index
    /// before `add_file` checks for `<file>.tbi` (:156).
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_tabix_fixture(
            dir.path(),
            "no_info.vcf.bgz",
            &["##fileformat=VCFv4.2", "#CHROM\tPOS", "chrM\t100"],
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("no_info.vcf.bgz.tbi")).unwrap();
        let mut plugin = GnomadMtPlugin::new();
        let message = expect_init_error(plugin.init(&[format!("file={}", file.display())]));
        assert!(
            message.contains("no_info.vcf.bgz.tbi not found"),
            "{message}"
        );
    }

    /// The `##INFO` lines are read as the bytes tabix prints them (gnomADMt.pm:146):
    /// a byte that is not UTF-8 (Latin-1 `é`, 0xE9) in a description leaves the
    /// field readable, the byte written as U+FFFD in the description.
    #[test]
    fn info_lines_that_are_not_utf8_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_bgzf_bytes(
            dir.path(),
            "gnomad.vcf.bgz",
            b"##fileformat=VCFv4.2\n##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"Fr\xE9quence\">\n\
              ##INFO=<ID=AF_het,Number=1,Type=Float,Description=\"Allele frequency\">\n#CHROM\tPOS\n",
        );
        assert_eq!(
            info_fields(&file).unwrap(),
            [
                ("AF_hom".to_string(), "Fr\u{FFFD}quence".to_string()),
                ("AF_het".to_string(), "Allele frequency".to_string()),
            ]
        );
    }

    #[test]
    fn file_without_af_hom_or_af_het_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"x\">",
            "##INFO=<ID=AC_hom,Number=1,Type=Integer,Description=\"y\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chrM\t3\t.\tT\tC\t.\tPASS\tAF_hom=0.5;AC_hom=1",
        ];
        let file = write_tabix_fixture(dir.path(), "no_het.vcf.bgz", &lines, IndexSpec::standard());
        let mut plugin = GnomadMtPlugin::new();
        assert_eq!(
            expect_init_error(plugin.init(&[format!("file={}", file.display())])),
            format!(
                "gnomADMt: {} does not carry the AF_hom and AF_het INFO fields; give the file as gnomAD distributes it",
                file.display()
            )
        );

        let mut plugin = GnomadMtPlugin::new();
        assert_eq!(
            expect_init_error(plugin.init(&[])),
            "gnomADMt: no file specified; give file=path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
        );
    }

    /// `gnomADMt.pm:146` opens the file through ensembl-io's VCF parser, which
    /// ends the run on a `##fileformat` line naming a major version other than
    /// 4 (`BaseVCF4.pm:125`) or no `<digits>.<digits>` version (`:129`), on
    /// whichever header line carries it; a `VCFv4.5` line and a file with no
    /// such line are read (`_parse_metadata_line` runs on a `##key=value` line
    /// alone).
    #[test]
    fn a_fileformat_line_outside_vcf_4_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let info = [
            "##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"x\">",
            "##INFO=<ID=AF_het,Number=1,Type=Float,Description=\"y\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chrM\t3\t.\tT\tC\t.\tPASS\tAF_hom=0.5;AF_het=0.1",
        ];
        let init = |name: &str, first: Option<&str>| {
            let lines: Vec<&str> = first.into_iter().chain(info).collect();
            let file = write_tabix_fixture(dir.path(), name, &lines, IndexSpec::standard());
            let mut plugin = GnomadMtPlugin::new();
            let result = plugin.init(&[
                format!("file={}", file.display()),
                "fields=AF_hom".to_string(),
            ]);
            (file, result)
        };
        let (file, result) = init("v33.vcf.bgz", Some("##fileformat=VCFv3.3"));
        assert_eq!(
            expect_init_error(result),
            format!(
                "gnomADMt: {}: ##fileformat=VCFv3.3: The VCF file format version 3.3 is not \
                 compatible with the parser version (VCF v4.2)",
                file.display()
            )
        );
        let (file, result) = init("v4.vcf.bgz", Some("##fileformat=VCFv4"));
        assert_eq!(
            expect_init_error(result),
            format!(
                "gnomADMt: {}: ##fileformat=VCFv4: The script can't read the VCF file format \
                 version of 'fileformat'",
                file.display()
            )
        );
        assert!(init("v45.vcf.bgz", Some("##fileformat=VCFv4.5")).1.is_ok());
        assert!(init("none.vcf.bgz", None).1.is_ok());
    }

    /// The descriptions are the file's `##INFO` descriptions as written, the
    /// comparison signs included.
    #[test]
    fn default_fields_are_the_five_haplogroup_fields() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        assert_eq!(
            keys(&plugin),
            [
                "gnomAD_hap_AC_hom",
                "gnomAD_hap_AC_het",
                "gnomAD_hap_AF_hom",
                "gnomAD_hap_AF_het",
                "gnomAD_hap_AN"
            ]
        );
        assert_eq!(
            description(&plugin, "gnomAD_hap_AN"),
            format!("List of overall allele number for each haplogroup, {HAP_ORDER}")
        );
        let plugin = init_with(dir.path(), &["fields=AF_het:AC_hom"]).unwrap();
        assert_eq!(keys(&plugin), ["gnomAD_AF_het", "gnomAD_AC_hom"]);
        assert_eq!(
            description(&plugin, "gnomAD_AF_het"),
            "Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < 0.95"
        );
    }

    #[test]
    fn fields_all_expands_to_every_info_id() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=all"]).unwrap();
        assert_eq!(
            keys(&plugin),
            [
                "gnomAD_filters",
                "gnomAD_hap_defining_variant",
                "gnomAD_AC_hom",
                "gnomAD_AF_hom",
                "gnomAD_AF_het",
                "gnomAD_hap_AN",
                "gnomAD_hap_AC_het",
                "gnomAD_hap_AC_hom",
                "gnomAD_hap_AF_hom",
                "gnomAD_hap_AF_het",
                "gnomAD_hap_hl_hist",
                "gnomAD_pop_AN",
                "gnomAD_pop_AC_hom",
                "gnomAD_pop_hl_hist",
            ]
        );
        // A field the matched row lacks, or carries as a flag, is written with
        // the empty value Ensembl VEP writes for an undefined one.
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_filters").as_deref(), Some(""));
        assert_eq!(
            value(&variants[0], "gnomAD_hap_defining_variant").as_deref(),
            Some("")
        );
        assert_eq!(value(&variants[0], "gnomAD_AC_hom").as_deref(), Some("19"));
        let variants = annotate(&plugin, vec![variant("MT", 22, 22, "T", "-")]);
        assert_eq!(
            value(&variants[0], "gnomAD_filters").as_deref(),
            Some("npg")
        );
    }

    #[test]
    fn unknown_fields_warn_and_all_unknown_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom:nope:AF_hom"]).unwrap();
        assert_eq!(keys(&plugin), ["gnomAD_AC_hom", "gnomAD_AF_hom"]);
        let message = init_message(dir.path(), &["fields=nope:zilch"]);
        assert!(
            message.starts_with(
                "gnomADMt: all fields given are invalid. Available fields are: filters, hap_defining_variant, AC_hom"
            ),
            "{message}"
        );
    }

    /// `hap_filter=HV:A` keeps the HV and A entries of every haplogroup field,
    /// in that order, and states that order in their descriptions (`:120-130`,
    /// `:226-233`, `:340-351`).
    #[test]
    fn hap_filter_reorders_values_and_rewrites_the_description() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["hap_filter=HV:A"]).unwrap();
        assert_eq!(
            description(&plugin, "gnomAD_hap_AN"),
            "List of overall allele number for each haplogroup, haplogroup order: ['HV','A']"
        );
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AN").as_deref(),
            Some("30|10")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AF_het").as_deref(),
            Some("0.3|0.1")
        );

        // An unknown haplogroup is skipped with a warning; none known is an error.
        let plugin = init_with(dir.path(), &["hap_filter=B:ZZ"]).unwrap();
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_hap_AN").as_deref(), Some("20"));
        let message = init_message(dir.path(), &["hap_filter=ZZ"]);
        assert_eq!(
            message,
            "gnomADMt: all haplogroups given are invalid. Available haplogroups are: A, B, HV"
        );

        // Without a haplogroup field selected the filter is warned about and ignored.
        let plugin = init_with(dir.path(), &["fields=AC_hom", "hap_filter=HV"]).unwrap();
        assert!(plugin.haplogroups.is_none());
    }

    /// Ensembl VEP 116.2 reads `hap_filter=all` and `pop_filter=all` as every
    /// group the file's order list states, in file order, before any name is
    /// validated (`gnomADMt.pm:82`, called from `:174` and `:195`), so the
    /// values keep their order and the descriptions are rewritten with the
    /// whole list in the filter's own spelling (`:126`).
    #[test]
    fn filter_all_keeps_every_group_in_file_order_and_rewrites_the_description() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["hap_filter=all"]).unwrap();
        assert_eq!(
            description(&plugin, "gnomAD_hap_AN"),
            "List of overall allele number for each haplogroup, haplogroup order: ['A','B','HV']"
        );
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AN").as_deref(),
            Some("10|20|30")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AF_het").as_deref(),
            Some("0.1|0.2|0.3")
        );

        let plugin = init_with(
            dir.path(),
            &["fields=pop_AC_hom:pop_hl_hist:hap_AN", "pop_filter=all"],
        )
        .unwrap();
        assert_eq!(
            description(&plugin, "gnomAD_pop_AC_hom"),
            "List of AC_hom for each population, population order: ['afr','eas','fin']"
        );
        assert_eq!(
            description(&plugin, "gnomAD_hap_AN"),
            format!("List of overall allele number for each haplogroup, {HAP_ORDER}")
        );
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(
            value(&variants[0], "gnomAD_pop_AC_hom").as_deref(),
            Some("7|8|9")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_pop_hl_hist").as_deref(),
            Some("1|1|1,2|2|2,3|3|3")
        );
    }

    #[test]
    fn pop_filter_reorders_values_and_rewrites_the_description() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(
            dir.path(),
            &["fields=pop_AC_hom:pop_AN:hap_AN", "pop_filter=fin:afr"],
        )
        .unwrap();
        assert_eq!(
            description(&plugin, "gnomAD_pop_AN"),
            "List of overall allele number for each population, population order: ['fin','afr']"
        );
        assert_eq!(
            description(&plugin, "gnomAD_hap_AN"),
            format!("List of overall allele number for each haplogroup, {HAP_ORDER}")
        );
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(
            value(&variants[0], "gnomAD_pop_AN").as_deref(),
            Some("300|100")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_pop_AC_hom").as_deref(),
            Some("9|7")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AN").as_deref(),
            Some("10|20|30")
        );
    }

    /// The `*_hl_hist` fields list one histogram per group, `,`-joined, so a
    /// filter re-indexes them on the comma and every other group field on `|`
    /// (`:344`, `:357`).
    #[test]
    fn hl_hist_fields_split_on_comma() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(
            dir.path(),
            &[
                "fields=hap_hl_hist:pop_hl_hist:hap_AN",
                "hap_filter=HV:B",
                "pop_filter=eas",
            ],
        )
        .unwrap();
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(
            value(&variants[0], "gnomAD_hap_hl_hist").as_deref(),
            Some("0|0|3,0|2|0")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_pop_hl_hist").as_deref(),
            Some("2|2|2")
        );
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AN").as_deref(),
            Some("30|20")
        );
        assert_eq!(
            description(&plugin, "gnomAD_hap_hl_hist"),
            "Histogram of heteroplasmy levels for each haplogroup; bin edges are: [0.0, 0.5, 1.0], haplogroup order: ['HV','B']"
        );
    }

    /// The file names the contig `chrM`; a record on `MT`, `M` or `chrM` reads
    /// it (`:256-261`).
    #[test]
    fn mt_and_chrm_contigs_both_match() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom"]).unwrap();
        for chr in ["MT", "M", "chrM", "chrMT"] {
            let variants = annotate(&plugin, vec![variant(chr, 3, 3, "T", "C")]);
            assert_eq!(
                value(&variants[0], "gnomAD_AC_hom").as_deref(),
                Some("19"),
                "{chr}"
            );
        }
        let variants = annotate(&plugin, vec![variant("21", 3, 3, "T", "C")]);
        assert!(variants[0].plugin_data.is_empty());
    }

    /// A deletion `AT>A` at 21 and an insertion `T>TC` at 42 are written with
    /// the shared base before them; the record's own trimmed alleles match
    /// them (`:306-314`), and an allele the file lacks matches nothing.
    #[test]
    fn shared_first_base_indels_are_trimmed() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom"]).unwrap();
        let variants = annotate(
            &plugin,
            vec![
                variant("MT", 22, 22, "T", "-"),
                variant("MT", 43, 42, "-", "C"),
                variant("MT", 43, 42, "-", "G"),
                variant("MT", 8993, 8993, "T", "G"),
                variant("MT", 8993, 8993, "T", "C"),
                variant("MT", 8993, 8993, "T", "A"),
            ],
        );
        let ac: Vec<Option<String>> = variants.iter().map(|v| value(v, "gnomAD_AC_hom")).collect();
        assert_eq!(
            ac,
            [
                Some("2".to_string()),
                Some("88".to_string()),
                None,
                Some("0".to_string()),
                Some("1".to_string()),
                None
            ]
        );
    }

    /// Ensembl VEP 116.2 trims the record's alleles as well as the row's inside
    /// `get_matched_variant_alleles` (ensembl-variation Utils/Sequence.pm:1199-1208
    /// and :1223-1232), so the substitution at 3 written with its unchanged
    /// neighbour on either side (`AT/AC` at 2, `TC/CC` at 3) takes the `T C`
    /// row at 3; written with another base changing, it takes nothing.
    #[test]
    fn a_substitution_with_an_unchanged_flanking_base_matches_its_single_base_row() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom"]).unwrap();
        let variants = annotate(
            &plugin,
            vec![
                variant("MT", 2, 3, "AT", "AC"),
                variant("MT", 3, 4, "TC", "CC"),
                variant("MT", 2, 3, "AT", "AG"),
            ],
        );
        let ac: Vec<Option<String>> = variants.iter().map(|v| value(v, "gnomAD_AC_hom")).collect();
        assert_eq!(ac, [Some("19".to_string()), Some("19".to_string()), None]);
    }

    /// Ensembl VEP 116.2 strips a row's shared first base, when its alleles
    /// differ in length, before `get_matched_variant_alleles` reads the row
    /// (gnomADMt.pm:306-314), so the row `AA A` at 500 is `A/-` at 501, which
    /// the deletion at 501 takes and the deletion at 500 does not, and the row
    /// `AA AAA` at 600 is `A/AA` at 601, which the insertions at 601 and 602
    /// take and the insertion at 600 does not.
    #[test]
    fn a_row_is_stripped_of_its_shared_first_base_before_the_match() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom"]).unwrap();
        let variants = annotate(
            &plugin,
            vec![
                variant("MT", 501, 501, "A", "-"),
                variant("MT", 500, 500, "A", "-"),
                variant("MT", 601, 600, "-", "A"),
                variant("MT", 602, 601, "-", "A"),
                variant("MT", 600, 599, "-", "A"),
            ],
        );
        let ac: Vec<Option<String>> = variants.iter().map(|v| value(v, "gnomAD_AC_hom")).collect();
        assert_eq!(
            ac,
            [
                Some("7".to_string()),
                None,
                Some("5".to_string()),
                Some("5".to_string()),
                None
            ]
        );
    }

    /// The values are written at record level, so a record with no transcript
    /// row carries them, and a record's transcript rows carry none of their own.
    #[test]
    fn intergenic_variants_are_annotated() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        let mut with_rows = variant("MT", 3, 3, "T", "C");
        with_rows
            .transcript_consequences
            .push(TranscriptConsequence {
                feature_type: FeatureType::Transcript,
                ..TranscriptConsequence::default()
            });
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C"), with_rows]);
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AC_hom").as_deref(),
            Some("4|5|6")
        );
        assert_eq!(
            value(&variants[1], "gnomAD_hap_AC_hom").as_deref(),
            Some("4|5|6")
        );
        assert!(variants[1].transcript_consequences[0]
            .plugin_data
            .is_empty());
    }

    #[test]
    fn hap_filter_without_an_order_list_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"x\">",
            "##INFO=<ID=AF_het,Number=1,Type=Float,Description=\"y\">",
            "##INFO=<ID=hap_AN,Number=1,Type=String,Description=\"List of overall allele number for each haplogroup\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chrM\t3\t.\tT\tC\t.\tPASS\tAF_hom=0.5;AF_het=0.1;hap_AN=1|2",
        ];
        let file = write_tabix_fixture(
            dir.path(),
            "no_order.vcf.bgz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = GnomadMtPlugin::new();
        let result = plugin.init(&[
            format!("file={}", file.display()),
            "fields=hap_AN".to_string(),
            "hap_filter=A".to_string(),
        ]);
        assert_eq!(
            expect_init_error(result),
            format!(
                "gnomADMt: hap_filter needs the haplogroup order list in the hap_AN INFO description of {}",
                file.display()
            )
        );
    }

    /// Ensembl VEP splits the line htslib hands it, the line end alone removed,
    /// on tabs (gnomADMt.pm:304) and the INFO column's pairs on `=`
    /// (`:327-329`), so the whitespace before the line end stays on the value
    /// of the column's last pair: `hap_AN=10|20|30 ` prints `10|20|30 ` and
    /// `hap_AN= ` prints ` `.
    #[test]
    fn a_trailing_space_on_the_last_info_pair_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF_hom,Number=1,Type=Float,Description=\"x\">",
            "##INFO=<ID=AF_het,Number=1,Type=Float,Description=\"y\">",
            "##INFO=<ID=hap_AN,Number=1,Type=String,Description=\"List of overall allele number for each haplogroup, haplogroup order: ['A', 'B', 'HV']\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chrM\t70\t.\tA\tG\t.\tPASS\tAF_hom=0.5;AF_het=0.1;hap_AN=10|20|30 ",
            "chrM\t80\t.\tA\tG\t.\tPASS\tAF_hom=0.5;AF_het=0.1;hap_AN= ",
        ];
        let file = write_tabix_fixture(
            dir.path(),
            "trailing_space.vcf.bgz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = GnomadMtPlugin::new();
        plugin
            .init(&[
                format!("file={}", file.display()),
                "fields=hap_AN".to_string(),
            ])
            .unwrap();
        let variants = annotate(
            &plugin,
            vec![
                variant("MT", 70, 70, "A", "G"),
                variant("MT", 80, 80, "A", "G"),
            ],
        );
        assert_eq!(
            value(&variants[0], "gnomAD_hap_AN").as_deref(),
            Some("10|20|30 ")
        );
        assert_eq!(value(&variants[1], "gnomAD_hap_AN").as_deref(), Some(" "));
    }

    #[test]
    fn header_pairs_keep_quoted_commas_and_signs() {
        let pairs = header_pairs(
            "ID=x,Number=1,Type=String,Description=\"a, b >= c \\\"q\\\" [d]\",Source=\"s\"",
        );
        let expected: Vec<(&str, String)> = vec![
            ("ID", "x".to_string()),
            ("Number", "1".to_string()),
            ("Type", "String".to_string()),
            ("Description", "a, b >= c \"q\" [d]".to_string()),
            ("Source", "s".to_string()),
        ];
        assert_eq!(pairs, expected);
        assert_eq!(
            group_order("x, haplogroup order: ['A', 'L0', 'HV'] tail", "haplogroup"),
            Some(vec!["A".to_string(), "L0".to_string(), "HV".to_string()])
        );
        assert_eq!(group_order("no list here", "haplogroup"), None);
        assert!(is_group_field("hap_AN", "hap_"));
        assert!(is_group_field("faf_hap_x", "hap_"));
        assert!(!is_group_field("hap_", "hap_"));
        assert!(!is_group_field("hapmax_AF_hom", "hap_"));
    }

    /// Ensembl VEP 116.2 never hands gnomADMt a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2220). A `<DEL>` over the file's
    /// records gets no field, and nor does a matching substitution carrying the
    /// structural mark, which pins the gate rather than the allele match; the
    /// plain substitution is annotated.
    #[test]
    fn a_structural_variant_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["fields=AC_hom"]).unwrap();
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = variant("chrM", 3, 50, "N", "<DEL>");
        deletion.is_structural = true;
        let mut flagged = variant("chrM", 3, 3, "T", "C");
        flagged.is_structural = true;
        let variants = annotate(
            &plugin,
            vec![deletion, flagged, variant("chrM", 3, 3, "T", "C")],
        );
        assert!(variants[0].plugin_data.is_empty());
        assert!(variants[1].plugin_data.is_empty());
        assert_eq!(value(&variants[2], "gnomAD_AC_hom").as_deref(), Some("19"));
    }

    /// Ensembl VEP puts every row of the window to `get_matched_variant_alleles`
    /// in file order until one matches (gnomADMt.pm:275-296), and the match
    /// dies on a row whose reference or alternate allele is undefined
    /// (`trim_sequences`, ensembl-variation Utils/Sequence.pm:968, reached from
    /// `:1230`): the plugin's `went wrong` warning names the row and the
    /// variant's line carries none of the fields, while a variant whose window
    /// holds no such row keeps its fields. `:304` assigns `split /\t/` to a
    /// list of eight scalars, which keeps trailing empty cells, so the alternate
    /// is undefined only on a row of four cells, the reference on one of three.
    #[test]
    fn an_undefined_allele_row_before_the_match_fails_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(dir.path(), "short_before.vcf.bgz", &["chrM\t3\t.\tT"], true);
        let plugin = init_file(&file);
        let variants = annotate_failing_on_an_undefined_row(
            &plugin,
            vec![
                variant("MT", 3, 3, "T", "C"),
                variant("MT", 8993, 8993, "T", "C"),
            ],
        );
        assert!(
            variants[0].plugin_data.is_empty(),
            "{:?}",
            variants[0].plugin_data
        );
        assert_eq!(value(&variants[1], "gnomAD_AC_hom").as_deref(), Some("1"));

        let file = fixture_with_rows(dir.path(), "no_ref_before.vcf.bgz", &["chrM\t3\t."], true);
        let plugin = init_file(&file);
        let variants =
            annotate_failing_on_an_undefined_row(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert!(variants[0].plugin_data.is_empty());
    }

    /// The loop reads on past a row that does not match (`:292-296`), so a row
    /// with no alternate cell behind the row for another allele at 3 fails the
    /// allele that row does not match, and the window opens two bases before
    /// the start for an insertion as for any variant (`:265-267`), so such a
    /// row at 41 fails the insertion at 43 whose own row sits at 42.
    #[test]
    fn an_undefined_allele_row_after_a_row_that_does_not_match_fails_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(dir.path(), "short_after.vcf.bgz", &["chrM\t3\t.\tT"], false);
        let plugin = init_file(&file);
        let variants =
            annotate_failing_on_an_undefined_row(&plugin, vec![variant("MT", 3, 3, "T", "G")]);
        assert!(variants[0].plugin_data.is_empty());

        let file = fixture_with_rows(dir.path(), "short_at_41.vcf.bgz", &["chrM\t41\t.\tT"], true);
        let plugin = init_file(&file);
        let variants =
            annotate_failing_on_an_undefined_row(&plugin, vec![variant("MT", 43, 42, "-", "C")]);
        assert!(variants[0].plugin_data.is_empty());
    }

    /// The loop stops at its first match (`last`, `:294`), so a row with no
    /// alternate cell behind the matching row is never read and the variant
    /// keeps its fields.
    #[test]
    fn an_undefined_allele_row_after_the_match_is_never_reached() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(dir.path(), "short_after.vcf.bgz", &["chrM\t3\t.\tT"], false);
        let plugin = init_file(&file);
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_AC_hom").as_deref(), Some("19"));
    }

    /// The window is `[start - 2, end]`, the end raised to the start for an
    /// insertion (`:265-267`): a row with no alternate cell at 40 lies outside
    /// the window of the insertion at 43 and is never read.
    #[test]
    fn an_undefined_allele_row_outside_the_window_is_never_reached() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(dir.path(), "short_at_40.vcf.bgz", &["chrM\t40\t.\tT"], true);
        let plugin = init_file(&file);
        let variants = annotate(&plugin, vec![variant("MT", 43, 42, "-", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_AC_hom").as_deref(), Some("88"));
    }

    /// The VCF preset of tabix ends a row at the last base of its reference
    /// allele, so a row at 40 whose reference spans three bases is returned for
    /// the window `[41, 43]` of the insertion at 43 and, lacking its alternate
    /// cell, dies there, its own `INFO` cell being absent changing nothing of
    /// that span; a two-base reference at 39 ends at 40, outside the window.
    #[test]
    fn a_short_row_whose_reference_reaches_the_window_fails_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(dir.path(), "span_into.vcf.bgz", &["chrM\t40\t.\tTAC"], true);
        let plugin = init_file(&file);
        let variants =
            annotate_failing_on_an_undefined_row(&plugin, vec![variant("MT", 43, 42, "-", "C")]);
        assert!(variants[0].plugin_data.is_empty());

        let file = fixture_with_rows(dir.path(), "span_short.vcf.bgz", &["chrM\t39\t.\tTA"], true);
        let plugin = init_file(&file);
        let variants = annotate(&plugin, vec![variant("MT", 43, 42, "-", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_AC_hom").as_deref(), Some("88"));
    }

    /// `:304` assigns `split /\t/` to a list of eight scalars, a split with a
    /// limit that keeps trailing empty cells, so an empty alternate cell is the
    /// defined empty string whether a cell follows it or not: the row dies
    /// nowhere, trims to `-`, matches no substitution, and the row behind it
    /// supplies the fields.
    #[test]
    fn an_empty_alt_cell_is_a_defined_row_that_matches_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let file = fixture_with_rows(
            dir.path(),
            "empty_alt_before.vcf.bgz",
            &[
                "chrM\t3\t.\tT\t",
                "chrM\t3\t.\tT\t\t.\tPASS\tAC_hom=0;AF_hom=0.5;AF_het=0.1",
            ],
            true,
        );
        let plugin = init_file(&file);
        let variants = annotate(&plugin, vec![variant("MT", 3, 3, "T", "C")]);
        assert_eq!(value(&variants[0], "gnomAD_AC_hom").as_deref(), Some("19"));
    }

    /// A row's POS cell reaches the match as `_parse_vcf` leaves it
    /// (`gnomADMt.pm:306-314`): untouched for a substitution, so the key holds
    /// its text (ensembl-variation Utils/Sequence.pm:1232) and `1555.0` and
    /// `1555.5` annotate nothing where `1555` does, nor does `1556x` at 1556;
    /// stepped by `$start++` for an indel sharing its first base, so `21.0` and
    /// `6x` are 22 and 7, the variants' starts, while `42.5` is 43.5 and
    /// matches nothing.
    #[test]
    fn a_position_cell_is_matched_as_parse_vcf_and_the_key_read_it() {
        let dir = tempfile::tempdir().unwrap();
        let ac = |n: u32| format!("AC_hom={n};AF_hom=0.5;AF_het=0.25");
        let mut lines: Vec<String> = rows()
            .into_iter()
            .filter(|line| line.starts_with('#'))
            .collect();
        lines.extend([
            format!("chrM\t6x\t.\tC\tCCTCAA\t.\tPASS\t{}", ac(61)),
            format!("chrM\t21.0\t.\tAT\tA\t.\tPASS\t{}", ac(62)),
            format!("chrM\t42.5\t.\tT\tTC\t.\tPASS\t{}", ac(63)),
            format!("chrM\t1555\t.\tA\tG\t.\tPASS\t{}", ac(64)),
            format!("chrM\t1555.0\t.\tA\tC\t.\tPASS\t{}", ac(65)),
            format!("chrM\t1555.5\t.\tA\tT\t.\tPASS\t{}", ac(66)),
            format!("chrM\t1556x\t.\tA\tG\t.\tPASS\t{}", ac(67)),
        ]);
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let file = write_tabix_fixture(dir.path(), "cells.vcf.bgz", &refs, IndexSpec::standard());
        let plugin = init_file(&file);
        let ac_hom = |start: u64, end: u64, reference: &str, alt: &str| {
            let variants = annotate(&plugin, vec![variant("MT", start, end, reference, alt)]);
            value(&variants[0], "gnomAD_AC_hom")
        };
        assert_eq!(ac_hom(7, 6, "-", "CTCAA").as_deref(), Some("61"));
        assert_eq!(ac_hom(22, 22, "T", "-").as_deref(), Some("62"));
        assert_eq!(ac_hom(43, 42, "-", "C"), None);
        assert_eq!(ac_hom(1555, 1555, "A", "G").as_deref(), Some("64"));
        assert_eq!(ac_hom(1555, 1555, "A", "C"), None);
        assert_eq!(ac_hom(1555, 1555, "A", "T"), None);
        assert_eq!(ac_hom(1556, 1556, "A", "G"), None);

        // A row short of its allele cells fails the variant, unless its position
        // cell is Perl-false, which returns the match first (Sequence.pm:1173).
        let short = |pos: &str| TabixRecord {
            chr: "chrM".to_string(),
            start: 1555,
            end: 1555,
            columns: vec!["chrM".into(), pos.into(), ".".into()],
            header: None,
        };
        let snv = variant("MT", 1555, 1555, "A", "G");
        let matching = |pos: &str| plugin.matching_values(&[&short(pos)], &snv);
        assert!(matches!(matching("1555"), Err(PluginError::Run(_))));
        assert!(matches!(matching("0"), Ok(None)));
        assert!(matches!(matching(""), Ok(None)));
    }
}
