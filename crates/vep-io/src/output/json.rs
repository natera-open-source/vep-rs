// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! JSON output format: one JSON object per input record, one line each (JSON Lines).
//!
//! The object is VEP's: the record's identity and coordinates at the top level,
//! then one `<feature type>_consequences` array per feature type whose entries
//! are the same field hash the text formats print, with keys lowercased, VEP's
//! renames applied (`Consequence` to `consequence_terms`, `Feature` to
//! `transcript_id`, `Allele` to `variant_allele`, ...), position strings split
//! into `*_start`/`*_end`, `YES` to `1`, and numeric-looking values emitted as
//! numbers except for the identifier keys.

use std::io::Write;

use indexmap::IndexMap;
use serde_json::Value;

use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::custom_annotation::{
    CustomAnnotations, CustomFieldValue, CustomRecord, CustomStats, CustomValue,
};
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;
use vep_core::variation::perl_number_string;

use crate::error::IoError;

use super::fields::{self, rendered_terms, FieldOptions, DEFAULT_OUTPUT_COLS};
use super::frequency::{self, FrequencyGroups};
use super::OutputFormatter;

/// JSON output formatter.
pub struct JsonOutputFormatter {
    options: FieldOptions,
    /// [`fields::vf_level_flag_fields`] of `options`: the keys written on the
    /// record object and left off every consequence entry.
    vf_level_fields: Vec<&'static str>,
    plugin_fields: Vec<String>,
    assembly: Option<String>,
}

impl JsonOutputFormatter {
    /// `assembly` is the name Ensembl VEP writes as `assembly_name`; `None`
    /// omits the key. The `vep` binary passes the resolved cache's: the
    /// cache's own, else the one its directory name carries, else the
    /// `--assembly` flag's.
    pub fn new(
        mut options: FieldOptions,
        plugin_fields: Vec<String>,
        assembly: Option<String>,
    ) -> Self {
        options.set_plugin_fields(plugin_fields.iter().cloned());
        let vf_level_fields = fields::vf_level_flag_fields(&options);
        Self {
            options,
            vf_level_fields,
            plugin_fields,
            assembly,
        }
    }

    /// The JSON value for one input record (all alleles split from one line).
    /// Every object, an unannotated record's included, carries `vep_rs_version`;
    /// a record with no allele is the object carrying that key alone.
    pub fn format_record_value(&self, record: &[&InputVariant]) -> Value {
        let mut root: IndexMap<String, Value> = IndexMap::new();
        root.insert(
            "vep_rs_version".to_string(),
            Value::String(fields::VEP_RS_VERSION.to_string()),
        );
        let Some(&first) = record.first() else {
            return Value::Object(root.into_iter().collect());
        };
        if record.iter().any(|v| v.vep_skip) {
            // VEP emits the unannotated record as its input line alone.
            if let Some(raw) = &first.raw_input {
                root.insert("input".to_string(), Value::String(raw.to_string()));
            }
            return Value::Object(root.into_iter().collect());
        }
        root.insert("id".to_string(), Value::String(first.variation_name()));
        root.insert(
            "seq_region_name".to_string(),
            Value::String(first.original_chr.clone()),
        );
        let start = first.original_start.unwrap_or(first.start);
        let end = first.original_end.unwrap_or(first.end);
        root.insert("start".to_string(), Value::Number(start.into()));
        root.insert("end".to_string(), Value::Number(end.into()));
        root.insert(
            "strand".to_string(),
            Value::Number((first.strand.as_i8() as i64).into()),
        );
        // A breakend record lists its own breakend first; the record string is
        // the paired form's.
        let string_source = record
            .iter()
            .find(|v| !v.is_single_breakend)
            .copied()
            .unwrap_or(first);
        let allele_string = string_source.working_allele_string();
        root.insert(
            "allele_string".to_string(),
            Value::String(allele_string.clone()),
        );
        if let Some(asm) = &self.assembly {
            root.insert("assembly_name".to_string(), Value::String(asm.clone()));
        }
        if let Some(raw) = &first.raw_input {
            root.insert("input".to_string(), Value::String(raw.to_string()));
        }
        if let Some(custom) = custom_annotations_object(record) {
            root.insert("custom_annotations".to_string(), custom);
        }
        // The keys VEP fills from the variant alone sit on the record object,
        // lower-cased, and on no consequence entry (`OutputFactory/JSON.pm`
        // copies `VariationFeature_to_output_hash` onto the record and builds
        // the entries from an empty hash).
        for name in &self.vf_level_fields {
            if *name == "NEAREST" {
                continue;
            }
            let key = json_key(name);
            // The structural variants stay a list (`[sort keys ...]`,
            // OutputFactory.pm:962-964), each name numberified (JSON.pm:288), and
            // so does the `--individual_zyg` list of `sample:ZYG` entries
            // (OutputFactory.pm:981-987); the `:` in each entry keeps it a string
            // under numberify.
            let value = match (*name, &first.sample_zygosities) {
                ("SV", _) => {
                    if first.overlapping_svs.is_empty() {
                        continue;
                    }
                    numberified_list(&first.overlapping_svs)
                }
                ("ZYG", Some(samples)) => {
                    if samples.entries.is_empty() {
                        continue;
                    }
                    Value::Array(
                        samples
                            .entries
                            .iter()
                            .map(|sample| Value::String(sample.to_string()))
                            .collect(),
                    )
                }
                // The record's hash holds no key of a `--custom` source
                // (`OutputFactory/JSON.pm:249-258`), so a source's values under
                // a key that spells the field are no value of the field.
                _ if fields::value_comes_from_a_custom_source(name, first, None, &self.options) => {
                    continue;
                }
                _ => {
                    let value = fields::field_value_without_custom_pushes(
                        name,
                        first,
                        None,
                        &self.options,
                        ",",
                    );
                    if value.is_empty() {
                        continue;
                    }
                    numberify(&key, value)
                }
            };
            root.insert(key, value);
        }
        // `nearest` is a record key that is a list (`JSON.pm:262` writes the
        // record hash's array as is) and is kept when empty, so the loop above
        // skips it and it is written here as an array, the empty list included.
        if let Some(nearest) = &first.nearest {
            root.insert(
                "nearest".to_string(),
                Value::Array(nearest.iter().map(|s| Value::String(s.clone())).collect()),
            );
        }

        let mut buckets: IndexMap<String, Vec<Value>> = IndexMap::new();
        let mut all_terms: Vec<Consequence> = Vec::new();
        for (variant, tc) in fields::record_rows(record) {
            // A summary row with no term is the record's hash alone, written here
            // as the object with no consequence entry; Ensembl VEP's JSON factory
            // dereferences the hash's absent consequence list and dies on it
            // (OutputFactory/JSON.pm:374), the divergence
            // `docs/divergences/json-summary-without-term.md` lists.
            if tc.is_none() && variant.summary_terms.as_ref().is_some_and(|t| t.is_empty()) {
                continue;
            }
            let bucket = match tc {
                Some(tc) => {
                    all_terms.extend(tc.consequences.iter().copied());
                    format!("{}_consequences", feature_type_key(tc))
                }
                None => {
                    all_terms.push(if variant.intergenic_overlap {
                        Consequence::IntergenicVariant
                    } else {
                        variant
                            .most_severe_consequence
                            .unwrap_or(Consequence::IntergenicVariant)
                    });
                    "intergenic_consequences".to_string()
                }
            };
            let obj = self.consequence_object(variant, tc);
            buckets.entry(bucket).or_default().push(Value::Object(obj));
        }
        for (bucket, entries) in buckets {
            root.insert(bucket, Value::Array(entries));
        }
        let most_severe = Consequence::most_severe(&all_terms)
            .map(|c| self.options.terms.term(c).to_string())
            .unwrap_or_else(|| "?".to_string());
        root.insert(
            "most_severe_consequence".to_string(),
            Value::String(most_severe),
        );

        if !first.colocated_variants.is_empty() {
            let groups = self.options.frequency_groups().for_json();
            let cvs: Vec<Value> = first
                .colocated_variants
                .iter()
                .map(|cv| Value::Object(colocated_variant_object(cv, &allele_string, &groups)))
                .collect();
            root.insert("colocated_variants".to_string(), Value::Array(cvs));
        }
        Value::Object(serde_json::Map::from_iter(root))
    }

    /// The JSON value for a single allele treated as a whole record. Embedding
    /// consumers that annotate one allele at a time use this.
    pub fn format_variant_value(&self, variant: &InputVariant) -> Value {
        self.format_record_value(&[variant])
    }

    fn consequence_object(
        &self,
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
    ) -> serde_json::Map<String, Value> {
        if let (None, Some(terms)) = (tc, &variant.summary_terms) {
            // `summary_only` is given an empty hash in the JSON factory
            // (`add_VariationFeatureOverlapAllele_info`, JSON.pm:320), so the
            // record's one entry carries its consequence terms alone; the
            // variant-level keys sit at the top level.
            let consequence_terms = rendered_terms(terms, self.options.terms, true)
                .into_iter()
                .map(|term| Value::String(term.to_string()))
                .collect();
            return serde_json::Map::from_iter([(
                "consequence_terms".to_string(),
                Value::Array(consequence_terms),
            )]);
        }
        let mut obj: IndexMap<String, Value> = IndexMap::new();
        // Under `--shift_genomic 1` VEP writes `Location` into every transcript
        // row's hash of a sequence variant (OutputFactory.pm:1711-1715; a
        // structural variant's rows take `:2127-2134`, which writes none) and
        // into an intergenic row's when `--shift_3prime 1` gave its insertion or
        // deletion a shift length, 0 included (:2000-2007), so those entries
        // carry a `location`.
        let include_location = match tc {
            Some(row) => {
                self.options.shift_genomic
                    && !variant.is_structural
                    && row.feature_type == vep_core::consequence::FeatureType::Transcript
            }
            None => {
                self.options.shift_3prime
                    && self.options.shift_genomic
                    && variant.intergenic_shift.is_some()
            }
        };
        let mut names: Vec<&str> = DEFAULT_OUTPUT_COLS
            .iter()
            .copied()
            .filter(|c| match *c {
                "Uploaded_variation" | "Existing_variation" => false,
                "Location" => include_location,
                _ => true,
            })
            .collect();
        // A plugin key named as a column the object leaves out prints under its
        // own name on the rows the plugin wrote it, as every key of the
        // consequence hash does (OutputFactory/JSON.pm:321-343; the keys left
        // out are the record's, `%SKIP_KEYS` at :93-96 and :262).
        if fields::plugin_names_a_fixed_column(variant, tc) {
            for c in DEFAULT_OUTPUT_COLS {
                if !names.contains(&c)
                    && fields::plugin_value(c, variant, tc, &self.options).is_some()
                {
                    names.push(c);
                }
            }
        }
        // The co-located variant and frequency fields are stubbed out of the
        // consequence hash by the JSON factory (OutputFactory/JSON.pm:493-515) and
        // written under `colocated_variants` instead; `SV` is the record's key
        // (OutputFactory.pm:962-964), written at the top level as
        // `colocated_structural_variants` and in no consequence object.
        // A plugin that groups its records for this format writes them as one
        // object per record under its name, and its keys leave the flat fields.
        let has_groups = !variant.plugin_records.is_empty();
        let grouped = |key: &str| {
            has_groups
                && variant
                    .plugin_records
                    .values()
                    .flatten()
                    .any(|record| record.contains_key(key))
        };
        // The JSON factory moves the custom annotations to the record before
        // any row hash is built (OutputFactory/JSON.pm:249-258), so a row's
        // `source` is its source cache alone (OutputFactory.pm:1592): a
        // `--custom` source named `SOURCE` lists its names in the text formats
        // and in `custom_annotations` here, on no row. The same holds for a
        // source's key that spells any other flag field: a row carries the
        // field's own value or none.
        let custom_source_alone = |key: &str| {
            (key == "SOURCE"
                && tc.is_none_or(|row| row.source_cache().is_none())
                && fields::plugin_value(key, variant, tc, &self.options).is_none())
                || fields::value_comes_from_a_custom_source(key, variant, tc, &self.options)
        };
        for k in fields::extra_keys(variant, tc, &self.options) {
            if !names.contains(&k)
                && !self.vf_level_fields.contains(&k)
                && !fields::colocated_field(k)
                && !grouped(k)
                && !custom_source_alone(k)
            {
                names.push(k);
            }
        }
        // A declared field named as one of the thirteen columns is already
        // listed, or left out with the column: the core value never stands in
        // for a plugin's unwritten key on a consequence object. Nor does a
        // `--custom` source's value under a key the plugin declares
        // (OutputFactory/JSON.pm:249-258), nor the record's or a co-located
        // variant's value under a declared field named as a field of the
        // record's hash, a co-located variant field or a frequency column,
        // which no row's hash holds
        // (OutputFactory/JSON.pm:261-262, 320, 493-515).
        for p in &self.plugin_fields {
            if !DEFAULT_OUTPUT_COLS.contains(&p.as_str())
                && !names.contains(&p.as_str())
                && !grouped(p)
                && !custom_source_alone(p)
                && !fields::value_is_held_by_no_json_row(p, variant, tc, &self.options)
            {
                names.push(p.as_str());
            }
        }
        for (group, records) in &variant.plugin_records {
            obj.insert(
                json_key(group),
                Value::Array(records.iter().map(plugin_record_object).collect()),
            );
        }
        // The JSON factory lower-cases every key of the row's hash
        // (OutputFactory/JSON.pm:345) and writes the keys it derives over the
        // result: the two ends of a position (:333-343), each prediction
        // tool's parts (:352-360) and the renamed keys (:377-383). They are
        // held here and written after the row's own names, so a plugin column
        // whose lower-cased name is one of them gives way; `None` is an end
        // the split leaves out.
        let mut derived: Vec<(String, Option<Value>)> = Vec::with_capacity(16);
        for name in names {
            // The SpliceVault fields are one nested object (SpliceVault.pm:191),
            // built from the whole row when the first of them the row carries is
            // reached, from the values as the plugin wrote them (the events a
            // list, never the text formats' `,` join, so no text form is
            // rendered for them); inside it an undefined cell is `null` and an
            // empty one `""`, since `OutputFactory/JSON.pm:327` reads the
            // top-level keys alone and the encoder writes the nested hash as it
            // stands.
            if SPLICEVAULT_FIELDS.contains(&name) {
                if !obj.contains_key(SPLICEVAULT_KEY)
                    && fields::plugin_value(name, variant, tc, &self.options).is_some()
                {
                    let object = splicevault_object(|field| {
                        match fields::plugin_value(field, variant, tc, &self.options) {
                            Some(text) if !text.is_empty() => Some(text.to_string()),
                            _ if fields::empty_plugin_value_is_a_string(field, variant, tc) => {
                                Some(String::new())
                            }
                            _ => None,
                        }
                    });
                    obj.insert(SPLICEVAULT_KEY.to_string(), object);
                }
                continue;
            }
            // The SpliceAI parts and flag are one nested object (SpliceAI.pm:316),
            // built from the whole row when the first of them the row carries is
            // reached, and the AlphaMissense columns another
            // (AlphaMissense.pm:281-282), one key at a time; inside either an
            // undefined value is `null` and an empty one `""`, as the nested hash
            // is encoded as it stands.
            if SPLICEAI_FIELDS.contains(&name) {
                if !obj.contains_key(SPLICEAI_KEY)
                    && fields::plugin_value(name, variant, tc, &self.options).is_some()
                {
                    let object = spliceai_object(|field| {
                        nested_plugin_value(field, variant, tc, &self.options)
                    });
                    obj.insert(SPLICEAI_KEY.to_string(), object);
                }
                continue;
            }
            if name.starts_with(ALPHAMISSENSE_FIELD_PREFIX) {
                if let Some(value) = nested_plugin_value(name, variant, tc, &self.options) {
                    let group = obj
                        .entry(ALPHAMISSENSE_KEY.to_string())
                        .or_insert_with(|| Value::Object(serde_json::Map::new()));
                    if let Value::Object(columns) = group {
                        columns.insert(name.to_string(), numberify_nested_value(name, value));
                    }
                }
                continue;
            }
            let value =
                fields::field_value_without_custom_pushes(name, variant, tc, &self.options, ",");
            // The two lists a regulatory or motif row carries empty stay empty
            // arrays (`OutputFactory/JSON.pm:327` drops an undefined value alone).
            if matches!(name, "TRANSCRIPTION_FACTORS" | "CELL_TYPE")
                && value.is_empty()
                && fields::empty_value_is_written(name, variant, tc, &self.options)
            {
                obj.insert(json_key(name), Value::Array(Vec::new()));
                continue;
            }
            // A plugin's empty string is a value (`OutputFactory/JSON.pm:327`
            // drops an undefined value alone).
            if value.is_empty()
                && fields::empty_plugin_value_is_a_string(name, variant, tc)
                && fields::empty_plugin_value_is_written(name, variant, tc, &self.options)
            {
                write_key(
                    &mut obj,
                    &mut derived,
                    name.to_lowercase(),
                    Value::String(String::new()),
                );
                continue;
            }
            // So is the empty reference of a feature string that opens with `/`.
            if value.is_empty()
                && name == "REF_ALLELE"
                && fields::empty_value_is_written(name, variant, tc, &self.options)
            {
                obj.insert(json_key(name), Value::String(String::new()));
                continue;
            }
            if value.is_empty() || (value == "-" && !fields::dash_is_a_value(name)) {
                continue;
            }
            match name {
                "Consequence" => {
                    let terms = value
                        .split(',')
                        .map(|t| Value::String(t.to_string()))
                        .collect();
                    derived.push(("consequence_terms".to_string(), Some(Value::Array(terms))));
                }
                "Feature_type" | "NEAREST" => {}
                "Feature" => {
                    let key = format!(
                        "{}_id",
                        tc.map(feature_type_key)
                            .unwrap_or_else(|| "intergenic".to_string())
                    );
                    derived.push((key, Some(Value::String(value))));
                }
                "cDNA_position" | "CDS_position" | "Protein_position" => {
                    let coord = name.trim_end_matches("_position").to_lowercase();
                    let (start, end) = position_bounds(&value);
                    derived.push((format!("{coord}_start"), start));
                    derived.push((format!("{coord}_end"), end));
                }
                "SIFT" | "PolyPhen" => {
                    let tool = name.to_lowercase();
                    let (prediction, score) = parse_prediction_score(&value);
                    if let Some(prediction) = prediction {
                        derived.push((
                            format!("{tool}_prediction"),
                            Some(Value::String(prediction.to_string())),
                        ));
                    }
                    if let Some(score) = score {
                        let key = format!("{tool}_score");
                        let number = numberify(&key, score.to_string());
                        derived.push((key, Some(number)));
                    }
                }
                "DOMAINS" => {
                    let domains =
                        domains_json(tc.map(|x| x.domains.as_slice()).unwrap_or_default());
                    obj.insert("domains".to_string(), Value::Array(domains));
                }
                // The fields VEP's output hash holds as lists (OutputFactory.pm:1487,
                // 1549, 1557, 1581, 1601-1609, 1671 at release/116.3), encoded as
                // arrays; the one undefined element VEP can produce, the miRNA
                // structure read one past its end, is `null`, while an empty element
                // of a split accession list stays an empty string.
                "FLAGS"
                | "MANE"
                | "RefSeq"
                | "REFSEQ_MATCH"
                | "SWISSPROT"
                | "TREMBL"
                | "UNIPARC"
                | "UNIPROT_ISOFORM"
                | "miRNA"
                | "TRANSCRIPTION_FACTORS" => {
                    let empty_element_is_null = name == "miRNA";
                    let elements = value
                        .split(',')
                        .map(|element| {
                            if element.is_empty() && empty_element_is_null {
                                Value::Null
                            } else {
                                Value::String(element.to_string())
                            }
                        })
                        .collect();
                    write_key(
                        &mut obj,
                        &mut derived,
                        name.to_lowercase(),
                        Value::Array(elements),
                    );
                }
                // The list `get_cell_types` builds (OutputFactory.pm:1968-1978),
                // one `<name>:<activity>` element per entry: an epigenome name
                // may hold a comma, so the elements come from the row, never
                // from a split of the joined text.
                "CELL_TYPE" => {
                    let entries = tc
                        .and_then(|row| row.cell_type.as_deref())
                        .unwrap_or_default();
                    let elements = entries
                        .iter()
                        .map(|entry| Value::String(entry.to_string()))
                        .collect();
                    obj.insert(json_key(name), Value::Array(elements));
                }
                "GA4GH_VRS" => {
                    // The object is nested rather than encoded, and `numberify`
                    // makes both interval bounds numbers (OutputFactory.pm:1409,
                    // JSON.pm:288); the text formats keep `start` a string.
                    if let Some(vrs) = fields::row_ga4gh_vrs(variant, tc, &self.options) {
                        obj.insert("ga4gh_vrs".to_string(), vrs_object(vrs));
                    }
                }
                "5UTR_annotation" => {
                    obj.insert(json_key(name), utr_annotation_object(&value));
                }
                _ => {
                    let lowercased = name.to_lowercase();
                    let v = if value == "YES" {
                        Value::Number(1.into())
                    } else {
                        let key = renamed_key(&lowercased).unwrap_or(lowercased.as_str());
                        numberify(key, value)
                    };
                    write_key(&mut obj, &mut derived, lowercased, v);
                }
            }
        }
        for (key, value) in derived {
            match value {
                Some(value) => {
                    obj.insert(key, value);
                }
                None => {
                    obj.swap_remove(&key);
                }
            }
        }
        serde_json::Map::from_iter(obj)
    }
}

/// Writes a field's value under its lower-cased name, or holds it in
/// `derived` under the key the JSON factory renames that name to
/// (`renamed_key`).
fn write_key(
    obj: &mut IndexMap<String, Value>,
    derived: &mut Vec<(String, Option<Value>)>,
    lowercased: String,
    value: Value,
) {
    match renamed_key(&lowercased) {
        Some(renamed) => derived.push((renamed.to_string(), Some(value))),
        None => {
            obj.insert(lowercased, value);
        }
    }
}

/// The record's `custom_annotations` object (`OutputFactory/JSON.pm:249-258`):
/// under each `--custom` short name the statistics object first, then every
/// record of the input record, `exact` records with the allele they matched,
/// all alleles together. Nothing when no source listed a record; a source with
/// statistics alone then lists them (:256 adds them to every source with
/// statistics). Each record is its name, its allele, its fields and its score as
/// the source read them, every value numberified as `numberify` recurses into
/// the hash (Utils.pm:225-243); a `FILTER` list is an array (`File/VCF.pm:222`),
/// an undefined `PC` (an `exact` match computes none, `VCF.pm:420`) `null`, and
/// statistics a source asked for but never computed `null`.
fn custom_annotations_object(record: &[&InputVariant]) -> Option<Value> {
    // Each allele carries its own copy of every source's entry, so a later
    // allele adds only what the earlier ones lack; within one allele every
    // entry is a distinct source's and is kept, equal records or not: two
    // sources under one short name each push their records (File.pm:437).
    let mut entries: Vec<&CustomAnnotations> = Vec::new();
    for allele in record {
        let before_this_allele = entries.len();
        for entry in &allele.custom_annotations {
            if !entries[..before_this_allele].contains(&entry) {
                entries.push(entry);
            }
        }
    }
    if entries.iter().all(|entry| entry.records.is_empty()) {
        return None;
    }
    let mut object: IndexMap<String, Value> = IndexMap::new();
    for entry in &entries {
        if object.contains_key(&*entry.short_name) {
            continue;
        }
        let same_name: Vec<&CustomAnnotations> = entries
            .iter()
            .copied()
            .filter(|other| other.short_name == entry.short_name)
            .collect();
        let mut list: Vec<Value> = Vec::new();
        match fields::custom_stats(same_name.iter().copied()) {
            Some(CustomStats::Values(values)) => list.push(Value::Object(
                values
                    .iter()
                    .map(|(stat, value)| (stat.clone(), custom_value_json(stat, value)))
                    .collect(),
            )),
            Some(CustomStats::Undefined) => list.push(Value::Null),
            None => {}
        }
        list.extend(
            same_name
                .iter()
                .flat_map(|entry| &entry.records)
                .map(custom_record_object),
        );
        if !list.is_empty() {
            object.insert(entry.short_name.to_string(), Value::Array(list));
        }
    }
    Some(Value::Object(serde_json::Map::from_iter(object)))
}

/// One custom record as Ensembl VEP's hash holds it (`File/VCF.pm:231-260`,
/// `File.pm:459-466`): `name`, `allele` when the record matched one, `fields`
/// when it carries any and `score` when the source read one.
fn custom_record_object(record: &CustomRecord) -> Value {
    let mut object: IndexMap<String, Value> = IndexMap::new();
    object.insert("name".to_string(), numberify("name", record.name.clone()));
    if let Some(allele) = &record.allele {
        object.insert(
            "allele".to_string(),
            numberify("allele", String::from_utf8_lossy(allele).into_owned()),
        );
    }
    if !record.fields.is_empty() {
        let fields = record.fields.iter().map(|(field, value)| {
            let value = match value {
                CustomFieldValue::Text(text) if field == "PC" && text.is_empty() => Value::Null,
                CustomFieldValue::Text(text) => numberify(field, text.clone()),
                CustomFieldValue::List(items) => numberified_list(items),
                CustomFieldValue::Undefined => Value::Null,
            };
            (field.clone(), value)
        });
        object.insert(
            "fields".to_string(),
            Value::Object(serde_json::Map::from_iter(fields)),
        );
    }
    if let Some(score) = &record.score {
        object.insert("score".to_string(), custom_value_json("score", score));
    }
    Value::Object(serde_json::Map::from_iter(object))
}

/// A score or statistic in the JSON format: the text numberified, a number as
/// Ensembl VEP writes it, an undefined value `null`.
fn custom_value_json(key: &str, value: &CustomValue) -> Value {
    match value {
        CustomValue::Text(text) => numberify(key, text.clone()),
        CustomValue::Number(number) => numberify(key, perl_number_string(*number)),
        CustomValue::Undefined => Value::Null,
    }
}

/// The key the AlphaMissense plugin's fields nest under in the JSON format: the
/// plugin returns `{AlphaMissense => {...}}` there (AlphaMissense.pm:281-282 at
/// VEP_plugins release/116), a consequence key the JSON factory lowers
/// (OutputFactory/JSON.pm:345). Every field the plugin writes carries the `am_`
/// prefix (AlphaMissense.pm:135, 152), under its own name inside the object.
const ALPHAMISSENSE_KEY: &str = "alphamissense";

const ALPHAMISSENSE_FIELD_PREFIX: &str = "am_";

/// One record a plugin groups under its name: the keys as the plugin wrote
/// them, each value numberified as `numberify` recurses into a nested hash
/// (Utils.pm:225-233 at release/116.3), `null` for a value the record lacks.
fn plugin_record_object(record: &vep_core::variant::PluginRecord) -> Value {
    Value::Object(
        record
            .iter()
            .map(|(key, value)| {
                let value = match value {
                    Some(value) => numberify(key, value.clone()),
                    None => Value::Null,
                };
                (key.clone(), value)
            })
            .collect(),
    )
}

/// The key the SpliceVault plugin's fields nest under in the JSON format: the
/// plugin returns `{SpliceVault => {...}}` there (SpliceVault.pm:191 at
/// VEP_plugins release/116), a consequence key the JSON factory lowers
/// (OutputFactory/JSON.pm:345).
const SPLICEVAULT_KEY: &str = "splicevault";

const SPLICEVAULT_TOP_EVENTS: &str = "SpliceVault_top_events";

/// The SpliceVault plugin's fields (SpliceVault.pm:142-148).
const SPLICEVAULT_FIELDS: [&str; 7] = [
    "SpliceVault_site_sample_count",
    "SpliceVault_site_max_depth",
    "SpliceVault_SpliceAI_delta",
    "SpliceVault_out_of_frame_events",
    "SpliceVault_site_pos",
    "SpliceVault_site_type",
    SPLICEVAULT_TOP_EVENTS,
];

/// The keys of one event's object, in the order `split(/:/, $event)` fills them
/// (SpliceVault.pm:184-188); the rank, the first part, keys the event.
const SPLICEVAULT_EVENT_KEYS: [&str; 4] = [
    "site_type",
    "transcript_impact",
    "percent_of_supporting_samples",
    "frame",
];

/// The SpliceVault plugin's fields as the JSON format nests them
/// (SpliceVault.pm:180-191): each field under its own name, the numeric-looking
/// values as numbers (`numberify` recurses into the plugin's hash,
/// OutputFactory/JSON.pm:288), and `SpliceVault_top_events` an object keyed by
/// each event's rank, the event's other parts under [`SPLICEVAULT_EVENT_KEYS`].
/// The events field carries the list of `:216`, one `rank:type:impact:percent:frame`
/// per element ([`fields::PLUGIN_LIST_SEPARATOR`]); a part the event lacks is
/// `null`, as the undefined value `split` leaves is, a trailing empty part is
/// `""`, and a part past the fifth is dropped. `field_text` gives a field's
/// text form on the row, `None` for a cell the table row lacks, which the
/// twelve-variable split of `SpliceVault.pm:201` leaves undefined and the
/// encoder writes as `null`; the events of such a cell are the empty list, `{}`
/// (`:216`).
fn splicevault_object(field_text: impl Fn(&str) -> Option<String>) -> Value {
    let mut object = serde_json::Map::new();
    for field in SPLICEVAULT_FIELDS {
        let text = field_text(field);
        let nested = if field == SPLICEVAULT_TOP_EVENTS {
            splicevault_events_object(text.as_deref().unwrap_or_default())
        } else {
            match text {
                Some(text) => numberify(field, text),
                None => Value::Null,
            }
        };
        object.insert(field.to_string(), nested);
    }
    Value::Object(object)
}

/// [`splicevault_object`]'s event objects keyed by rank, one per element of the
/// list the plugin carries, the `split(/\|/)` of the column (SpliceVault.pm:216)
/// with its `|` kept ([`fields::PLUGIN_LIST_SEPARATOR`]), so an event is one
/// whatever its text holds, a `,` included: an interior empty event (`a||b`) is
/// keyed `""` with every part `null`, as the undefined rank `split(/:/, '')`
/// leaves is, and the empty list is `{}`. The parts of an event are the
/// `split(/:/, $event)` of `:184` into the rank and the four keyed variables,
/// which runs with a limit one above that count and so keeps a trailing empty
/// part as the empty string where the list-context split of the column drops one.
fn splicevault_events_object(events: &str) -> Value {
    let variables = SPLICEVAULT_EVENT_KEYS.len() + 1;
    let mut by_rank = serde_json::Map::new();
    for event in vep_core::sv_routing::perl_split(events, &[fields::PLUGIN_LIST_SEPARATOR]) {
        let mut parts = event.splitn(variables + 1, ':');
        let rank = parts.next().unwrap_or_default().to_string();
        let mut parts_by_key = serde_json::Map::new();
        for key in SPLICEVAULT_EVENT_KEYS {
            let part = match parts.next() {
                Some(part) => numberify(key, part.to_string()),
                None => Value::Null,
            };
            parts_by_key.insert(key.to_string(), part);
        }
        by_rank.insert(rank, Value::Object(parts_by_key));
    }
    Value::Object(by_rank)
}

/// The key the SpliceAI plugin's fields nest under in the JSON format: the
/// plugin returns `{SpliceAI => {...}}` there (SpliceAI.pm:316 at VEP_plugins
/// release/116), a consequence key the JSON factory lowers
/// (OutputFactory/JSON.pm:345).
const SPLICEAI_KEY: &str = "spliceai";

/// The prefix the SpliceAI plugin's nine fields carry in the VCF format and
/// under `split_output`, and drop inside the JSON object (SpliceAI.pm:277-287).
const SPLICEAI_PREFIX: &str = "SpliceAI_pred_";

/// The SpliceAI plugin's fields of the nine-field shape (SpliceAI.pm:279-287)
/// and its cutoff flag (`:302`), which keeps its name inside the object.
const SPLICEAI_FIELDS: [&str; 10] = [
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
    "SpliceAI_cutoff",
];

/// The SpliceAI plugin's fields as the JSON format nests them (SpliceAI.pm:275-288,
/// 302, 316): each present field under its name without the `SpliceAI_pred_`
/// prefix, the numeric-looking values as numbers (`numberify` recurses into the
/// plugin's hash, OutputFactory/JSON.pm:288) except `SYMBOL`, which the factory
/// exempts (JSON.pm:120). A member that holds the values of several entries,
/// joined by commas, is no number and stays that text, as the member
/// `26,G` of Ensembl VEP's own object does. `nested` gives a field's value on
/// the row ([`nested_plugin_value`]); a field the row lacks is left out.
fn spliceai_object(nested: impl Fn(&str) -> Option<Value>) -> Value {
    let mut object = serde_json::Map::new();
    for field in SPLICEAI_FIELDS {
        let Some(value) = nested(field) else {
            continue;
        };
        let key = field.strip_prefix(SPLICEAI_PREFIX).unwrap_or(field);
        object.insert(key.to_string(), numberify_nested_value(key, value));
    }
    Value::Object(object)
}

/// A plugin field's value inside the object a plugin nests its fields in: its
/// text, for the caller to `numberify` under the nested key, the empty string
/// as `""` where the plugin wrote the empty string
/// ([`fields::empty_plugin_value_is_a_string`]) and as `null` where it wrote
/// its undefined value, since the encoder writes the nested hash as it stands
/// and `OutputFactory/JSON.pm:327` reads the top-level keys alone; `None` for a
/// field the row lacks.
fn nested_plugin_value(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> Option<Value> {
    let text = fields::plugin_value(field, variant, tc, options)?;
    Some(
        if !text.is_empty() || fields::empty_plugin_value_is_a_string(field, variant, tc) {
            Value::String(text.to_string())
        } else {
            Value::Null
        },
    )
}

/// A value of [`nested_plugin_value`] written under `key`: its text numberified
/// under that key, `null` as it stands.
fn numberify_nested_value(key: &str, value: Value) -> Value {
    match value {
        Value::String(text) => numberify(key, text),
        other => other,
    }
}

/// The UTRAnnotator plugin's `5UTR_annotation` as the JSON format nests it
/// (UTRAnnotator.pm:271-276 at VEP_plugins release/116). The text form, which
/// every text format prints, joins one `key=value` group per affected upstream
/// ORF with `,` and the pairs of a group with `:`; the object keys the groups
/// `1`, `2`, ... in that order, each a hash of its pairs with the
/// numeric-looking values as numbers (`numberify` recurses into the plugin's
/// hash).
fn utr_annotation_object(annotation: &str) -> Value {
    let groups = annotation.split(',').enumerate().map(|(index, group)| {
        let pairs = group.split(':').map(|pair| {
            let (key, value) = pair.split_once('=').unwrap_or((pair, ""));
            (key.to_string(), numberify(key, value.to_string()))
        });
        ((index + 1).to_string(), Value::Object(pairs.collect()))
    });
    Value::Object(groups.collect())
}

/// A GA4GH VRS allele as the JSON format nests it.
fn vrs_object(vrs: &vep_core::spdi::VrsAllele) -> Value {
    serde_json::json!({
        "location": {
            "type": "SequenceLocation",
            "interval": {"type": "SimpleInterval", "start": vrs.start, "end": vrs.end},
            "sequence_id": format!("refseq:{}", vrs.sequence_id),
        },
        "state": {"type": "SequenceState", "sequence": vrs.sequence},
        "type": "Allele",
    })
}

/// Splits a position column into the JSON `*_start` and `*_end` values as
/// `OutputFactory/JSON.pm` `add_VariationFeatureOverlapAllele_info` does: `N`
/// writes both, `N-M` one each, `N-?` the start alone, `?-N` the end alone. A
/// `--total_length` suffix (`N/T`) stays a string and is written at both ends.
/// A range keeps its end under the flag: `N-M/T` writes `N` and the string
/// `M/T`, `?-N/T` writes `N/T` under the end key alone and `N-?/T` the start
/// alone. Ensembl VEP writes `N` under both keys of `N-M/T` and neither key of
/// `?-N/T`: its test of a usable end (`JSON.pm:337`) admits no total, so the
/// end key takes the start. That is the documented difference
/// `total_length_json_range_end`.
fn position_bounds(value: &str) -> (Option<Value>, Option<Value>) {
    let (start, end) = value.split_once('-').unwrap_or((value, value));
    (position_value(start), position_value(end))
}

/// A position literal, `N` or `N/T`, as a JSON number or string; `None` for
/// anything else (`?`, empty).
fn position_value(literal: &str) -> Option<Value> {
    match literal.split_once('/') {
        None if is_digits(literal) => literal.parse::<u64>().ok().map(Value::from),
        Some((n, total)) if is_digits(n) && is_digits(total) => {
            Some(Value::String(literal.to_string()))
        }
        _ => None,
    }
}

fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn feature_type_key(tc: &TranscriptConsequence) -> String {
    // `RegulatoryFeature` -> `regulatory_feature`, `MotifFeature` -> `motif_feature`.
    tc.feature_type
        .as_str()
        .to_lowercase()
        .replace("feature", "_feature")
        .trim_start_matches('_')
        .to_string()
}

/// VEP's JSON key for a field name: lowercased, then renamed where VEP renames
/// the lowercased key (`renamed_key`), so a plugin's `chr` column is
/// `seq_region_name`.
fn json_key(field: &str) -> String {
    let key = field.to_lowercase();
    match renamed_key(&key) {
        Some(renamed) => renamed.to_string(),
        None => key,
    }
}

/// The key VEP renames a lowercased key to (`%RENAME_KEYS`,
/// OutputFactory/JSON.pm:98-112, :380-383), `None` for a key it keeps and for
/// `consequence` and `feature` (:99, :379), which `consequence_object` renames
/// in the arms of those two columns.
fn renamed_key(lowercased: &str) -> Option<&'static str> {
    Some(match lowercased {
        "gene" => "gene_id",
        "allele" => "variant_allele",
        "symbol" => "gene_symbol",
        "symbol_source" => "gene_symbol_source",
        "overlapbp" => "bp_overlap",
        "overlappc" => "percentage_overlap",
        "refseq" => "refseq_transcript_ids",
        "ensp" => "protein_id",
        "chr" => "seq_region_name",
        "variation_name" => "id",
        "sv" => "colocated_structural_variants",
        "clinical_impact" => "clinvar_somatic_classification",
        _ => return None,
    })
}

/// Numeric-looking values become numbers, as VEP's `numberify` does, except the
/// identifier keys it exempts (OutputFactory/JSON.pm:114-121; `SYMBOL` is the
/// SpliceAI plugin's gene symbol inside its nested object).
fn numberify(key: &str, value: String) -> Value {
    if matches!(
        key,
        "seq_region_name" | "id" | "gene_id" | "gene_symbol" | "transcript_id" | "SYMBOL"
    ) {
        return Value::String(value);
    }
    match perl_number(&value) {
        Some(number) => Value::Number(number),
        None => Value::String(value),
    }
}

/// The JSON number `numberify` (Utils.pm:231) writes for a value
/// `looks_like_number` admits, `None` for one it leaves as text. Blanks at
/// either end, the six characters of Perl's `isSPACE` (space, tab, line feed,
/// vertical tab, form feed, carriage return), are admitted and dropped by the
/// `+ 0`; between them an optionally signed integer, a decimal with digits on
/// one or both sides of the dot, or an exponent form in either case. An
/// integer is exact up to `u64::MAX`, and so is an exponent form naming a
/// whole number within the signed or unsigned 64-bit range, which Perl's
/// addition keeps in an integer slot; every other value is the double of its
/// 15-digit text, and a zero fraction below 1e15 writes as an integer
/// (`100.00` is 100), as `%.15g` prints a Perl number. The ten bytes
/// `0 but true`, Perl's idiom for a true zero, are the number 0 as they
/// stand and text with a blank at either end. A value whose double is not
/// finite (`Inf`, `NaN`, an exponent past the double's range) is `None`:
/// JSON has no such number, and Ensembl VEP writes the bare tokens `inf`
/// and `nan` there, which no JSON reader accepts.
fn perl_number(value: &str) -> Option<serde_json::Number> {
    let text =
        value.trim_matches(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r'));
    // Every finite number opens with a digit, a sign or a dot, so the three
    // parses are skipped for the identifiers and terms most values are; the
    // infinity and not-a-number words end as text either way, here when
    // unsigned and at the finite check below when signed.
    if !text
        .as_bytes()
        .first()
        .is_some_and(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.'))
    {
        return None;
    }
    if let Ok(i) = text.parse::<i64>() {
        return Some(i.into());
    }
    if let Ok(u) = text.parse::<u64>() {
        return Some(u.into());
    }
    let Some(f) = text.parse::<f64>().ok().filter(|f| f.is_finite()) else {
        // `looks_like_number` compares the whole value, blanks included, with
        // the ten bytes, so only the exact spelling is the number 0.
        return (value == "0 but true").then(|| 0.into());
    };
    if f.fract() == 0.0 && text.contains(['e', 'E']) {
        if (0.0..u64::MAX as f64).contains(&f) {
            return Some((f as u64).into());
        }
        if (i64::MIN as f64..0.0).contains(&f) {
            return Some((f as i64).into());
        }
    }
    // A decimal of at most 15 significant digits is recovered from its double
    // by a 15-digit print, so for such a text the format-and-parse would
    // return `f` unchanged and is skipped.
    let mantissa = text
        .split_once(['e', 'E'])
        .map_or(text, |(mantissa, _exponent)| mantissa);
    let significant_digits = mantissa
        .bytes()
        .filter(u8::is_ascii_digit)
        .skip_while(|b| *b == b'0')
        .count();
    let f = if significant_digits <= 15 {
        f
    } else {
        fifteen_significant_digits(f)
    };
    if f.fract() == 0.0 && f.abs() < 1e15 {
        return Some((f as i64).into());
    }
    serde_json::Number::from_f64(f)
}

/// The double of `f` written with 15 significant digits, the precision Perl
/// gives a numeric value's text (`%.15g`, which JSON::XS writes as well): a
/// value of 15 or fewer significant digits is itself, one with more is rounded
/// as Ensembl VEP's JSON output rounds it (`0.5500299334526062` writes
/// `0.550029933452606`). The four largest doubles round up past the double's
/// range, where Ensembl VEP's `1.79769313486232e+308` token reads back as an
/// infinity; they stay the finite value.
fn fifteen_significant_digits(f: f64) -> f64 {
    format!("{f:.14e}")
        .parse::<f64>()
        .ok()
        .filter(|rounded| rounded.is_finite())
        .unwrap_or(f)
}

/// Splits a `SIFT` or `PolyPhen` value into its prediction and score as
/// `OutputFactory/JSON.pm:353-360` (release/116.3) does with
/// `([a-z_]+)?\(?([\d.]+)?\)?`, matched case-insensitively: every part is
/// optional, so the match sits at the start of the value, the prediction is its
/// leading run of letters and underscores, and the score the run of digits and
/// dots after that run and an optional `(`. A bare score (`--sift s`) sets no
/// prediction and a bare prediction (`--sift p`) no score.
fn parse_prediction_score(value: &str) -> (Option<&str>, Option<&str>) {
    let letters = value
        .bytes()
        .take_while(|b| b.is_ascii_alphabetic() || *b == b'_')
        .count();
    let (prediction, rest) = value.split_at(letters);
    let rest = rest.strip_prefix('(').unwrap_or(rest);
    let digits = rest
        .bytes()
        .take_while(|b| b.is_ascii_digit() || *b == b'.')
        .count();
    (
        (!prediction.is_empty()).then_some(prediction),
        (digits > 0).then(|| &rest[..digits]),
    )
}

/// The `domains` array of a consequence as `OutputFactory/JSON.pm:363-371`
/// (release/116.3) builds it: each entry the `DOMAINS` column prints (the label
/// and identifier parts `grep {$_}` keeps, joined by `:`, OutputFactory.pm:1512)
/// re-split by [`word_colon_word`] into `{db, name}`. The split is Perl's
/// `m/(\w+)\:(\w+)/` followed by `push ... if $1 && $2`, and a failed match
/// leaves `$1` and `$2` at the values of the list's last successful one, so an
/// entry the pattern does not match repeats the previous entry's `db` and
/// `name`, and is dropped only when no earlier entry of the list matched. A
/// capture of `0` (`Gene3D:0.10.510.10` names `0`) is false to Perl, so that
/// entry pushes nothing, and neither does an unmatched entry after it until a
/// later match replaces the captures. `numberify` then runs over the whole
/// record (JSON.pm:288) and neither key is exempt, so a name of digits alone,
/// Gene3D's truncated `1`, is a number.
fn domains_json(domains: &[(String, String)]) -> Vec<Value> {
    let mut captured: Option<(String, String)> = None;
    let mut out = Vec::with_capacity(domains.len());
    for (label, id) in domains {
        let entry = match (perl_true(label), perl_true(id)) {
            (true, true) => format!("{label}:{id}"),
            (false, true) => id.clone(),
            (true, false) => label.clone(),
            (false, false) => String::new(),
        };
        if let Some(split) = word_colon_word(&entry) {
            captured = Some(split);
        }
        if let Some((db, name)) = captured
            .as_ref()
            .filter(|(db, name)| perl_true(db) && perl_true(name))
        {
            out.push(Value::Object(serde_json::Map::from_iter([
                ("db".to_string(), numberify("db", db.clone())),
                ("name".to_string(), numberify("name", name.clone())),
            ])));
        }
    }
    out
}

/// The `db` and `name` of a `DOMAINS` entry as `(\w+):(\w+)` reads them: the
/// first `:` with a word character on both sides splits the entry, the run of
/// word characters before it is the database and the run after it the name, so
/// `Gene3D:1.10.510.10` names `1` and `PANTHER:PTHR24416:SF8` names
/// `PTHR24416`; `None` for an entry with no such `:`. Word characters are ASCII
/// letters, digits and `_`.
fn word_colon_word(text: &str) -> Option<(String, String)> {
    let bytes = text.as_bytes();
    let is_word = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    let colon = (1..bytes.len().saturating_sub(1))
        .find(|&i| bytes[i] == b':' && is_word(bytes[i - 1]) && is_word(bytes[i + 1]))?;
    let db_start = bytes[..colon]
        .iter()
        .rposition(|b| !is_word(*b))
        .map_or(0, |i| i + 1);
    let name_end = bytes[colon + 1..]
        .iter()
        .position(|b| !is_word(*b))
        .map_or(bytes.len(), |i| colon + 1 + i);
    Some((
        text[db_start..colon].to_string(),
        text[colon + 1..name_end].to_string(),
    ))
}

/// The populations whose matched-allele frequencies a `colocated_variants`
/// entry nests under `frequencies` (`@FREQ_FIELDS`, OutputFactory/JSON.pm:130-136),
/// each written under its lower-cased name.
const FREQ_FIELDS: &[&str] = &[
    "AF",
    "AFR",
    "AMR",
    "ASN",
    "EAS",
    "SAS",
    "EUR",
    "AA",
    "EA",
    "ExAC",
    "ExAC_Adj",
    "ExAC_AFR",
    "ExAC_AMR",
    "ExAC_EAS",
    "ExAC_FIN",
    "ExAC_NFE",
    "ExAC_OTH",
    "ExAC_SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
];

/// One `colocated_variants` entry, as `add_colocated_variant_info_JSON`
/// (OutputFactory/JSON.pm:407-476 at release/116.3) writes the record's hash:
/// every cache column that is set, not empty and not `0` (:433-435), the
/// `failed` and matched-allele bookkeeping dropped, `chr` as `seq_region_name`
/// and `variation_name` as `id`
/// (:98-112), `clin_sig` and `pubmed` as the lists
/// [`ColocatedVariant::clin_sig_terms`](vep_core::variant::ColocatedVariant::clin_sig_terms)
/// and `pubmed_ids` split, `var_synonyms` as
/// `{source: [ids]}` split on `--`, `::` and `,` (:449-471) from the value
/// [`ColocatedVariant::synonyms`](vep_core::variant::ColocatedVariant::synonyms)
/// reads, `clinical_impact` as
/// the raw string under `clinvar_somatic_classification`, every other column the
/// cache declares (a pre-108 dump's `minor_allele` and `minor_allele_freq`)
/// under its own name, a `0` kept only under a name containing `af` (:415,
/// :433-435), and the frequency columns replaced by `frequencies`: for each
/// alternate allele of the record's allele string the frequencies the text
/// formats would print on that allele's rows under `groups` (JSON.pm:276-284
/// splits the record's `allele_string` and builds one frequency hash per allele
/// through `add_colocated_frequency_data`, so an allele of a `--minimal` split
/// record is read by its spelling in the record's string, not its trimmed one,
/// and matches nothing), `AF` under `af` and each other column under its
/// lower-cased name (:419-431). Every value is numberified except the
/// identifiers (:288, `%NUMBERIFY_EXEMPT`).
///
/// `allele_string` is the record's, as the object's parent prints it.
fn colocated_variant_object(
    cv: &vep_core::variant::ColocatedVariant,
    allele_string: &str,
    groups: &FrequencyGroups,
) -> serde_json::Map<String, Value> {
    let mut obj: IndexMap<String, Value> = IndexMap::new();
    if let Some(chr) = &cv.chr {
        obj.insert("seq_region_name".to_string(), Value::String(chr.clone()));
    }
    obj.insert("id".to_string(), Value::String(cv.id.clone()));
    obj.insert("start".to_string(), Value::Number(cv.start.into()));
    obj.insert("end".to_string(), Value::Number(cv.end.into()));
    obj.insert(
        "strand".to_string(),
        Value::Number((cv.strand as i64).into()),
    );
    if let Some(allele_string) = cv.allele_string.as_deref().filter(|s| perl_true(s)) {
        obj.insert(
            "allele_string".to_string(),
            numberify("allele_string", allele_string.to_string()),
        );
    }
    if cv.somatic {
        obj.insert("somatic".to_string(), Value::Number(1.into()));
    }
    if cv.phenotype {
        obj.insert("phenotype_or_disease".to_string(), Value::Number(1.into()));
    }
    if let Some(terms) = cv.clin_sig_terms() {
        obj.insert("clin_sig".to_string(), numberified_elements(terms));
    }
    if let Some(ids) = cv.pubmed_ids() {
        obj.insert("pubmed".to_string(), numberified_elements(ids));
    }
    if let Some(clinvar) = cv.clinvar.as_deref() {
        for (key, value) in [
            ("clin_sig_allele", &clinvar.clin_sig_allele),
            ("clin_sig_ref_allele", &clinvar.clin_sig_ref_allele),
            ("clinvar_somatic_classification", &clinvar.clinical_impact),
        ] {
            if let Some(value) = value.as_deref().filter(|v| perl_true(v)) {
                obj.insert(key.to_string(), numberify(key, value.to_string()));
            }
        }
    }
    if let Some(synonyms) = cv.synonyms().filter(|s| perl_true(s)) {
        obj.insert("var_synonyms".to_string(), var_synonyms_by_source(synonyms));
    }
    let mut other_columns: Vec<(&String, &String)> = cv
        .frequency_strings
        .iter()
        .filter(|(column, value)| {
            !FREQ_FIELDS.contains(&column.as_str())
                && (value.as_str() != "0" || column.contains("af"))
        })
        .collect();
    other_columns.sort();
    for (column, value) in other_columns {
        obj.insert(column.clone(), numberify(column, value.clone()));
    }
    let mut frequencies: IndexMap<String, Value> = IndexMap::new();
    for row_allele in allele_string.split('/').skip(1) {
        let mut by_population: IndexMap<String, Value> = IndexMap::new();
        for population in FREQ_FIELDS {
            if let Some(value) = frequency::reported_frequency(cv, population, row_allele, groups) {
                by_population.insert(
                    population.to_lowercase(),
                    numberify("frequency", value.text()),
                );
            }
        }
        if !by_population.is_empty() {
            frequencies.insert(
                row_allele.to_string(),
                Value::Object(serde_json::Map::from_iter(by_population)),
            );
        }
    }
    if !frequencies.is_empty() {
        obj.insert(
            "frequencies".to_string(),
            Value::Object(serde_json::Map::from_iter(frequencies)),
        );
    }
    serde_json::Map::from_iter(obj)
}

/// A `,`-split list field as the JSON factory writes it (OutputFactory/JSON.pm:469),
/// every element numberified (`Utils.pm:235-243` runs over array elements with no
/// exemption); a list of no element is written as `[]`.
fn numberified_list(values: &[String]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|value| numberify("", value.clone()))
            .collect(),
    )
}

/// [`numberified_list`] of the elements a stored list splits into.
fn numberified_elements<'a>(elements: impl Iterator<Item = &'a str>) -> Value {
    numberified_list(&elements.map(str::to_string).collect::<Vec<_>>())
}

/// `var_synonyms` as the JSON factory nests it (OutputFactory/JSON.pm:451-460):
/// the cache string split on `--` into `source::id,id` groups, each a key holding
/// the `,`-split ids of its second `::` field; a group without `::` holds no id,
/// and a repeated source keeps its last group.
fn var_synonyms_by_source(synonyms: &str) -> Value {
    let mut by_source: IndexMap<String, Value> = IndexMap::new();
    for group in perl_split(synonyms, "--") {
        let mut parts = perl_split(&group, "::").into_iter();
        let source = parts.next().unwrap_or_default();
        let ids = perl_split(&parts.next().unwrap_or_default(), ",");
        by_source.insert(source, numberified_list(&ids));
    }
    Value::Object(serde_json::Map::from_iter(by_source))
}

/// Perl's list-context `split` on a literal separator: trailing empty fields are
/// dropped, so an empty string yields no field.
fn perl_split(text: &str, separator: &str) -> Vec<String> {
    let mut fields: Vec<String> = text.split(separator).map(str::to_string).collect();
    while fields.last().is_some_and(String::is_empty) {
        fields.pop();
    }
    fields
}

impl OutputFormatter for JsonOutputFormatter {
    fn write_header(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        Ok(())
    }

    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError> {
        self.format_record(&[variant])
    }

    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        if record.is_empty() || record.iter().any(|v| v.oversize_sv) {
            return Ok(Vec::new());
        }
        let value = self.format_record_value(record);
        let line = serde_json::to_string(&value)
            .map_err(|e| IoError::VcfParse(format!("JSON serialization error: {}", e)))?;
        Ok(vec![line])
    }

    fn finish(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields::fixtures::{edit_transcript, with_rseq_codes};
    use serde_json::json;
    use smallvec::smallvec;
    use std::sync::Arc;
    use vep_core::consequence::{Consequence, FeatureType, Impact, TranscriptConsequence};
    use vep_core::terms::Vocabulary;
    use vep_core::transcript::Transcript;

    fn make_annotated_variant() -> InputVariant {
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.id = Some("rs123".into());
        variant.raw_input = Some("21\t25000100\trs123\tA\tG\t.\t.\t.".into());
        let mut tc = TranscriptConsequence::default();
        tc.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = "ENST00000352957".into();
            transcript.gene_stable_id = "ENSG00000141959".into();
            transcript.gene_symbol = Some("MRPL39".into());
            transcript.gene_symbol_source = Some("HGNC".into());
            transcript.hgnc_id = Some("HGNC:14027".into());
            transcript.biotype = "protein_coding".into();
            transcript.flags = Arc::from(["cds_end_NF".to_string()]);
            transcript
        }));
        tc.consequences = smallvec![Consequence::SynonymousVariant];
        tc.impact = Impact::LOW;
        tc.canonical = true;
        tc.cdna_position = Some("428".into());
        tc.cds_position = Some("366".into());
        tc.protein_position = Some("122".into());
        tc.amino_acids = Some("L".into());
        tc.codons = Some("ctG/ctA".into());
        tc.strand = -1;
        tc.feature_type = FeatureType::Transcript;
        variant.transcript_consequences = vec![tc];
        variant.most_severe_consequence = Some(Consequence::SynonymousVariant);
        variant
    }

    fn fmt(options: FieldOptions) -> JsonOutputFormatter {
        JsonOutputFormatter::new(options, Vec::new(), Some("GRCh37".into()))
    }

    /// A motif consequence as `OutputFactory/JSON.pm` writes it from the motif
    /// row's hash (`:311-391`): `high_inf_pos` stays the string `Y` or `N`
    /// (only `YES` becomes 1 and `numberify` leaves a letter alone, `:330`,
    /// `:288`), `motif_score_change` and `motif_pos` are numbers, `strand` a
    /// number, `transcription_factors` an array kept when empty, and the row
    /// sits under `motif_feature_consequences` with its `motif_feature_id`. A
    /// matrix without counts writes neither score key.
    #[test]
    fn motif_objects_match_json_pm() {
        let options = FieldOptions {
            regulatory: true,
            ..Default::default()
        };
        let motif_row = |name: &str,
                         high_inf_pos: Option<bool>,
                         score_change: Option<f64>,
                         factors: &[&str]| {
            let mut row = TranscriptConsequence::default();
            row.feature_id = Some(Arc::from("ENSM00000574814"));
            row.consequences = smallvec![Consequence::TfBindingSiteVariant];
            row.feature_type = FeatureType::MotifFeature;
            row.strand = -1;
            row.motif = Some(Box::new(vep_core::consequence::MotifFields::new(
                Arc::from(name),
                Some(15),
                high_inf_pos,
                score_change,
                factors.iter().map(|f| Arc::from(*f)).collect(),
            )));
            row
        };
        let mut variant = InputVariant::new(
            "21".into(),
            5_090_724,
            5_090_724,
            b"A".to_vec(),
            b"T".to_vec(),
        );
        variant.transcript_consequences = vec![
            motif_row("ENSPFM0484", Some(false), Some(0.0386), &["RUNX3"]),
            motif_row("ENSPFM0290", Some(true), Some(-0.0894), &[]),
            motif_row("ENSPFM_STUB", None, None, &["CTCF"]),
        ];
        variant.most_severe_consequence = Some(Consequence::TfBindingSiteVariant);
        let parsed = fmt(options).format_variant_value(&variant);
        let rows = parsed["motif_feature_consequences"].as_array().unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(
            rows[0],
            json!({
                "consequence_terms": ["TF_binding_site_variant"], "high_inf_pos": "N",
                "impact": "MODIFIER", "motif_feature_id": "ENSM00000574814",
                "motif_name": "ENSPFM0484", "motif_pos": 15, "motif_score_change": 0.039,
                "strand": -1, "transcription_factors": ["RUNX3"], "variant_allele": "T"
            })
        );
        assert_eq!(rows[1]["high_inf_pos"], "Y");
        assert_eq!(rows[1]["motif_score_change"], json!(-0.089));
        assert_eq!(rows[1]["transcription_factors"], json!([]));
        assert!(rows[2].get("high_inf_pos").is_none());
        assert!(rows[2].get("motif_score_change").is_none());
        assert_eq!(rows[2]["motif_name"], "ENSPFM_STUB");
        assert_eq!(parsed["most_severe_consequence"], "TF_binding_site_variant");
        assert!(parsed.get("transcript_consequences").is_none());
    }

    /// The reference allele of a feature string that opens with `/` is the empty
    /// string, which `OutputFactory/JSON.pm:327` keeps where it drops an
    /// undefined value: Ensembl VEP writes `"ref_allele":""` on the consequence
    /// of `21 30289183 x <empty> G` under `--show_ref_allele`.
    #[test]
    fn an_empty_reference_allele_is_written_as_an_empty_string() {
        let options = || FieldOptions {
            show_ref_allele: true,
            ..Default::default()
        };
        let mut v = make_annotated_variant();
        let parsed = fmt(options()).format_variant_value(&v);
        assert_eq!(parsed["transcript_consequences"][0]["ref_allele"], "A");
        v.allele_string = "/G".into();
        let parsed = fmt(options()).format_variant_value(&v);
        assert_eq!(parsed["transcript_consequences"][0]["ref_allele"], "");
    }

    #[test]
    fn top_level_matches_vep() {
        let v = make_annotated_variant();
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        assert_eq!(parsed["id"], "rs123");
        assert_eq!(parsed["seq_region_name"], "21");
        assert_eq!(parsed["start"], 25000100);
        assert_eq!(parsed["end"], 25000100);
        assert_eq!(parsed["strand"], 1);
        assert_eq!(parsed["allele_string"], "A/G");
        assert_eq!(parsed["assembly_name"], "GRCh37");
        assert_eq!(parsed["input"], "21\t25000100\trs123\tA\tG\t.\t.\t.");
        assert_eq!(parsed["most_severe_consequence"], "synonymous_variant");
        assert_eq!(parsed["vep_rs_version"], fields::VEP_RS_VERSION);
        assert!(
            parsed.get("variant_class").is_none(),
            "variant_class needs its flag"
        );
        assert!(parsed.get("existing_variation").is_none());
    }

    /// A numeric value carries at most 15 significant digits, the precision
    /// JSON::XS writes a Perl number with: a plugin score of 16 or 17 digits
    /// rounds to the double of its 15-digit text, a shorter value is itself, a
    /// value whose rounding is whole becomes an integer, and a frequency's
    /// exponent form stays the number it names.
    #[test]
    fn numbers_carry_fifteen_significant_digits_as_json_xs_writes_them() {
        assert_eq!(
            numberify("x", "0.5500299334526062".into()),
            json!(0.550029933452606)
        );
        assert_eq!(
            numberify("x", "0.27204954624176025".into()),
            json!(0.27204954624176)
        );
        assert_eq!(
            numberify("x", "0.9500000099999999".into()),
            json!(0.95000001)
        );
        assert_eq!(
            numberify("x", "0.84708172082901".into()),
            json!(0.84708172082901)
        );
        assert_eq!(numberify("x", "0.99999999999999999".into()), json!(1));
        assert_eq!(numberify("x", "6.875e-05".into()), json!(6.875e-5));
        assert_eq!(numberify("x", "-0.448027".into()), json!(-0.448027));
        assert_eq!(
            numberify("x", "12345678901234567".into()),
            json!(12345678901234567_i64)
        );
        assert_eq!(
            numberify("x", "1.7976931348623157e308".into()),
            json!(f64::MAX)
        );
        assert_eq!(
            numberify("x", format!("17976931348623157{}", "0".repeat(292))),
            json!(f64::MAX)
        );
    }

    /// `looks_like_number` (`Utils.pm:231`) admits blanks at either end of a
    /// value, the six ASCII characters Perl's `isSPACE` names, and `+ 0` drops
    /// them, so a padded cell is written as the number without its blanks; a
    /// blank between the sign and the digits, a value of blanks alone and a
    /// non-ASCII space keep the text a string.
    #[test]
    fn numberify_admits_the_blanks_looks_like_number_admits() {
        for padded in [
            "93 ", " 93", " 93 ", "\t93", "93\t", "93\r", "93\n", "93\x0b", "93\x0c", "  93  ",
        ] {
            assert_eq!(numberify("x", padded.into()), json!(93), "{padded:?}");
        }
        assert_eq!(numberify("x", " 0.63 ".into()), json!(0.63));
        assert_eq!(numberify("x", "-0.55399 ".into()), json!(-0.55399));
        assert_eq!(numberify("x", " 1e5".into()), json!(100000));
        assert_eq!(numberify("x", "+449 ".into()), json!(449));
        for text in [
            "+ 93", "- 93", " ", "  ", "\t", "", "93\u{a0}", "9 3", "1 000",
        ] {
            assert_eq!(numberify("x", text.into()), json!(text), "{text:?}");
        }
    }

    /// The shapes `looks_like_number` admits read as the number `+ 0` makes
    /// of them: a leading sign, a fraction without its integer or fraction
    /// digits, an exponent in either case with an optional sign, a leading
    /// zero as a decimal digit and a negative zero as 0; a lone sign or dot, a
    /// bare or trailing exponent marker, a dot after an exponent, hex, binary
    /// and octal prefixes, an underscore and a comma stay strings.
    #[test]
    fn numberify_reads_signs_fractions_and_exponents_as_looks_like_number_does() {
        let numbers = [
            ("+93", json!(93)),
            ("-93", json!(-93)),
            (".5", json!(0.5)),
            ("5.", json!(5)),
            ("+.5", json!(0.5)),
            ("-.5", json!(-0.5)),
            ("-5.", json!(-5)),
            ("1e5", json!(100000)),
            ("1E5", json!(100000)),
            ("1e+5", json!(100000)),
            ("1e-5", json!(1e-5)),
            ("1.5e3", json!(1500)),
            ("1.e5", json!(100000)),
            (".5e1", json!(5)),
            ("1e-400", json!(0)),
            ("017", json!(17)),
            ("00", json!(0)),
            ("007", json!(7)),
            ("-0", json!(0)),
            ("-0.0", json!(0)),
        ];
        for (text, number) in numbers {
            assert_eq!(numberify("x", text.into()), number, "{text:?}");
        }
        for text in [
            "+", "-", ".", "+.", "5..", ".5.", "1e", "e5", "1e5.0", "0x10", "0x1e", "0X10",
            "0b101", "0o17", "1_000", "1,000", "93x", "x93", "93%", "2-3", "chr21:1",
        ] {
            assert_eq!(numberify("x", text.into()), json!(text), "{text:?}");
        }
    }

    /// An integer past `i64` is still a JSON integer up to `u64::MAX`, as Perl
    /// keeps it exact in an unsigned value, and past that, or below `i64::MIN`,
    /// the double of its 15-digit text, the value Perl's numeric conversion
    /// gives it. An exponent form naming a whole number within either 64-bit
    /// range is that integer, exact, where a decimal form of the same value is
    /// the double `%.15g` prints: an integer below 1e15, an exponent form above.
    #[test]
    fn numberify_writes_an_integer_past_i64_as_perl_does() {
        let exponent_forms = [
            ("1e15", json!(1_000_000_000_000_000_i64)),
            ("1.5e15", json!(1_500_000_000_000_000_i64)),
            ("-1e15", json!(-1_000_000_000_000_000_i64)),
            ("1e16", json!(10_000_000_000_000_000_i64)),
            ("4.503599627370496e15", json!(4_503_599_627_370_496_i64)),
            (
                "1.2345678901234567e18",
                json!(1_234_567_890_123_456_768_i64),
            ),
            (
                "9.223372036854775807e18",
                json!(9_223_372_036_854_775_808_u64),
            ),
            ("2e64", json!(2e64)),
            ("-3e63", json!(-3e63)),
            ("2.9999999999999996e0", json!(3)),
            ("0e0", json!(0)),
            ("-0e0", json!(0)),
        ];
        for (text, number) in exponent_forms {
            assert_eq!(numberify("x", text.into()), number, "{text:?}");
        }
        assert_eq!(numberify("x", "2.9999999999999996".into()), json!(3));
        assert_eq!(
            numberify("x", "4503599627370496.0".into()),
            json!(4.5035996273705e15)
        );
        assert_eq!(numberify("x", "1000000000000000.0".into()), json!(1e15));
        assert_eq!(
            numberify("x", "10000000000000000".into()),
            json!(10_000_000_000_000_000_i64)
        );
        assert_eq!(
            numberify("x", "9223372036854775807".into()),
            json!(9223372036854775807_i64)
        );
        assert_eq!(
            numberify("x", "9223372036854775808".into()),
            json!(9223372036854775808_u64)
        );
        assert_eq!(
            numberify("x", "18446744073709551615".into()),
            json!(18446744073709551615_u64)
        );
        assert_eq!(
            numberify("x", "18446744073709551616".into()),
            json!(1.84467440737096e19)
        );
        assert_eq!(numberify("x", "99999999999999999999".into()), json!(1e20));
        assert_eq!(
            numberify("x", "-9223372036854775808".into()),
            json!(-9223372036854775808_i64)
        );
        assert_eq!(
            numberify("x", "-9223372036854775809".into()),
            json!(-9.22337203685478e18)
        );
    }

    /// `looks_like_number` admits the ten bytes `0 but true`, Perl's idiom
    /// for a true zero, exactly as they stand, and `+ 0` makes them 0; the
    /// words with a blank at either end, in another case, with a second
    /// blank or with another word stay the string.
    #[test]
    fn numberify_reads_0_but_true_as_the_number_0_as_looks_like_number_does() {
        assert_eq!(numberify("x", "0 but true".into()), json!(0));
        for text in [
            " 0 but true",
            "0 but true ",
            "0 but true\n",
            "0 But True",
            "0 but  true",
            "0 but false",
            "1 but true",
            "0 but true but",
        ] {
            assert_eq!(numberify("x", text.into()), json!(text), "{text:?}");
        }
    }

    /// A value whose number is not finite has no JSON number: infinity and
    /// not-a-number spellings, which `looks_like_number` admits, and an
    /// exponent past the double's range stay the string the text formats
    /// print. Ensembl VEP writes the bare tokens `inf` and `nan` for them,
    /// which are not JSON: the documented class `json_non_finite_number_token`.
    #[test]
    fn numberify_keeps_a_value_with_no_finite_number_as_text() {
        for text in [
            "Inf",
            "inf",
            "INF",
            "+Inf",
            "-Inf",
            "Infinity",
            "infinity",
            "-infinity",
            "NaN",
            "nan",
            "-nan",
            "NaNQ",
            "1.#INF",
            "1e400",
            "-1e400",
            " Inf ",
        ] {
            assert_eq!(numberify("x", text.into()), json!(text), "{text:?}");
        }
    }

    /// A plugin's grouped records are one array of objects per consequence
    /// entry under the plugin's lower-cased name, each object's keys as the
    /// plugin wrote them with numeric values numberified and a value the record
    /// lacks `null`; the group's `,`-joined flat fields, which the text formats
    /// print, are left out, and a plugin field outside any group stays flat.
    #[test]
    fn grouped_plugin_records_are_arrays_of_objects_in_place_of_their_flat_fields() {
        let mut v = make_annotated_variant();
        v.plugin_data
            .insert("OpenTargets_gwasGeneId".into(), "ENSG1,NA".into());
        v.plugin_data
            .insert("OpenTargets_gwasLocusToGeneScore".into(), "0.55,NA".into());
        // A made-up CADD_PHRED value.
        v.plugin_data.insert("CADD_PHRED".into(), "23.1".into());
        let mut gwas = vep_core::variant::PluginRecord::new();
        gwas.insert("OpenTargets_gwasGeneId".into(), Some("ENSG1".into()));
        gwas.insert(
            "OpenTargets_gwasLocusToGeneScore".into(),
            Some("0.55".into()),
        );
        let mut qtl = vep_core::variant::PluginRecord::new();
        qtl.insert("OpenTargets_gwasGeneId".into(), Some(String::new()));
        qtl.insert("OpenTargets_gwasLocusToGeneScore".into(), None);
        v.plugin_records
            .insert("OpenTargets".into(), vec![gwas, qtl]);
        let plugin_fields = vec![
            "CADD_PHRED".to_string(),
            "OpenTargets_gwasGeneId".to_string(),
            "OpenTargets_gwasLocusToGeneScore".to_string(),
        ];
        let formatter = JsonOutputFormatter::new(FieldOptions::default(), plugin_fields, None);
        let tc = &formatter.format_variant_value(&v)["transcript_consequences"][0];
        assert_eq!(
            tc["opentargets"],
            json!([
                {"OpenTargets_gwasGeneId": "ENSG1", "OpenTargets_gwasLocusToGeneScore": 0.55},
                {"OpenTargets_gwasGeneId": "", "OpenTargets_gwasLocusToGeneScore": null}
            ])
        );
        assert_eq!(tc["cadd_phred"], 23.1);
        assert!(tc.get("opentargets_gwasgeneid").is_none(), "{tc}");
        assert!(tc.get("opentargets_gwaslocustogenescore").is_none(), "{tc}");
    }

    /// A plugin's empty value stands for Ensembl VEP's undefined value, whose
    /// key the consequence object drops (OutputFactory/JSON.pm:327), except
    /// where the plugin defines an empty string: the GWAS plugin's empty
    /// columns, the Downstream plugin's `DownstreamProtein`, the translation
    /// cut at a stop in its first codon (Downstream.pm:139-147), and a table
    /// cell TranscriptAnnotator's row has and leaves empty, which the row names
    /// in `plugin_empty_strings` (TranscriptAnnotator.pm:245-250), keep their
    /// keys as the empty string; a row the plugin left without the key has no
    /// key, and so has a row whose table row lacks the cell.
    #[test]
    fn a_plugins_empty_value_is_kept_only_where_ensembl_vep_defines_it() {
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            vec![
                "DownstreamProtein".to_string(),
                "ProteinLengthChange".to_string(),
                "GWAS_pmid".to_string(),
                "Other".to_string(),
                "scores_SIFT_pred".to_string(),
                "scores_Comment".to_string(),
            ],
            None,
        );
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        tc.plugin_data
            .insert("DownstreamProtein".to_string(), String::new());
        tc.plugin_data
            .insert("ProteinLengthChange".to_string(), "-7".to_string());
        tc.plugin_data
            .insert("GWAS_pmid".to_string(), String::new());
        tc.plugin_data.insert("Other".to_string(), String::new());
        tc.plugin_data
            .insert("scores_SIFT_pred".to_string(), String::new());
        tc.plugin_data
            .insert("scores_Comment".to_string(), String::new());
        tc.plugin_empty_strings.push("scores_SIFT_pred".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["downstreamprotein"], json!(""));
        assert_eq!(row["proteinlengthchange"], -7);
        assert_eq!(row["gwas_pmid"], json!(""));
        assert!(row.get("other").is_none(), "{row}");
        assert_eq!(row["scores_sift_pred"], json!(""));
        assert!(row.get("scores_comment").is_none(), "{row}");
        let bare = formatter.format_variant_value(&make_annotated_variant())
            ["transcript_consequences"][0]
            .clone();
        assert!(bare.get("downstreamprotein").is_none(), "{bare}");
    }

    /// A plugin's keys are written over the line hash after the core fields
    /// (OutputFactory.pm:398, :2240), and the JSON factory renames the
    /// lowercased keys (`%RENAME_KEYS`, JSON.pm:98-112, :380-383): under
    /// `--plugin dbNSFP,<file>,ALL` the consequence object carries the file's
    /// `tsl` (a list in the made-up row, so a string) and its `chr` as `seq_region_name`, a string
    /// (`%NUMBERIFY_EXEMPT`, :114), on the rows the plugin wrote and no `tsl`
    /// on a row it left alone while `--tsl` is off.
    #[test]
    fn a_plugins_tsl_and_chr_print_as_ensembl_vep_names_them() {
        let plugin_fields = ["TSL", "chr", "SIFT_score"].map(str::to_string).to_vec();
        let formatter = JsonOutputFormatter::new(FieldOptions::default(), plugin_fields, None);
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        edit_transcript(tc, |t| t.tsl = Some("1".to_string()));
        tc.plugin_data
            .insert("TSL".to_string(), "5,1,2,3,2,NA".to_string());
        tc.plugin_data.insert("chr".to_string(), "21".to_string());
        tc.plugin_data
            .insert("SIFT_score".to_string(), "0.249".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["tsl"], json!("5,1,2,3,2,NA"));
        assert_eq!(row["seq_region_name"], json!("21"));
        assert!(row.get("chr").is_none(), "{row}");
        assert_eq!(row["sift_score"], json!(0.249));
        let mut bare = make_annotated_variant();
        edit_transcript(&mut bare.transcript_consequences[0], |t| {
            t.tsl = Some("3".to_string())
        });
        let bare = formatter.format_variant_value(&bare)["transcript_consequences"][0].clone();
        assert!(bare.get("tsl").is_none(), "{bare}");
        assert!(bare.get("seq_region_name").is_none(), "{bare}");
    }

    /// The JSON factory lower-cases every key of a row first
    /// (OutputFactory/JSON.pm:345), a plugin's `SIFT_score` to `sift_score`,
    /// and then writes each core tool's prediction and score over those keys,
    /// each only when the core value holds it (:352-360). So on a row with a
    /// core value the key carries the core's part, and the plugin's column
    /// stays where the core value has no such part: its prediction under
    /// `--sift s`, its score under `--sift p`, both on a row with no core
    /// value.
    #[test]
    fn a_core_prediction_score_and_a_plugin_column_of_the_same_key() {
        let columns = [
            ("SIFT_score", "0.249"),
            ("SIFT_prediction", "column_call"),
            ("PolyPhen_score", "0.5"),
            ("PolyPhen_prediction", "column_call"),
        ];
        let options = FieldOptions {
            sift: true,
            polyphen: true,
            ..Default::default()
        };
        let declared = columns.map(|(column, _)| column.to_string()).to_vec();
        let formatter = JsonOutputFormatter::new(options, declared, None);
        let row = |sift: Option<&str>, polyphen: Option<&str>| {
            let mut v = make_annotated_variant();
            let tc = &mut v.transcript_consequences[0];
            tc.sift = sift.map(str::to_string);
            tc.polyphen = polyphen.map(str::to_string);
            for (column, value) in columns {
                tc.plugin_data.insert(column.to_string(), value.to_string());
            }
            formatter.format_variant_value(&v)["transcript_consequences"][0].clone()
        };
        let keys = |row: &Value| {
            ["sift", "polyphen"].map(|tool| {
                (
                    row[format!("{tool}_prediction")].clone(),
                    row[format!("{tool}_score")].clone(),
                )
            })
        };
        // `--sift b --polyphen b`: the core value holds both parts.
        assert_eq!(
            keys(&row(Some("tolerated(1)"), Some("benign(0.09)"))),
            [
                (json!("tolerated"), json!(1)),
                (json!("benign"), json!(0.09))
            ]
        );
        // `s`: the core score, the column's prediction.
        assert_eq!(
            keys(&row(Some("0.01"), Some("0.998"))),
            [
                (json!("column_call"), json!(0.01)),
                (json!("column_call"), json!(0.998))
            ]
        );
        // `p`: the core prediction, the column's score.
        assert_eq!(
            keys(&row(Some("deleterious"), Some("probably_damaging"))),
            [
                (json!("deleterious"), json!(0.249)),
                (json!("probably_damaging"), json!(0.5))
            ]
        );
        // A row with one core value and a row with none.
        assert_eq!(
            keys(&row(Some("tolerated(1)"), None)),
            [
                (json!("tolerated"), json!(1)),
                (json!("column_call"), json!(0.5))
            ]
        );
        let no_core_value = row(None, None);
        assert_eq!(
            keys(&no_core_value),
            [
                (json!("column_call"), json!(0.249)),
                (json!("column_call"), json!(0.5))
            ]
        );
        assert!(no_core_value.get("sift").is_none(), "{no_core_value}");
    }

    /// The JSON factory writes the other keys it derives over the lower-cased
    /// ones too: the two ends of a split position (OutputFactory/JSON.pm:333-343)
    /// and the renamed keys, `consequence_terms`, `gene_id` and the feature's
    /// `<type>_id` among them (:377-383). A plugin column whose lower-cased
    /// name is such a key gives way on a row whose core field is split or
    /// renamed, an end the split leaves out included, and keeps the key on a
    /// row with no such core field.
    #[test]
    fn a_plugin_column_named_as_a_derived_key_gives_way_to_the_core_field() {
        let columns = [
            "consequence_terms",
            "Transcript_ID",
            "gene_id",
            "cdna_start",
            "cds_end",
            "protein_id",
            "gene_symbol",
        ];
        let declared = columns.map(str::to_string).to_vec();
        let formatter = JsonOutputFormatter::new(FieldOptions::default(), declared, None);
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        tc.cds_position = Some("366-?".into());
        for column in columns {
            tc.plugin_data
                .insert(column.to_string(), "column".to_string());
        }
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["consequence_terms"], json!(["synonymous_variant"]));
        assert_eq!(row["transcript_id"], json!("ENST00000352957"));
        assert_eq!(row["gene_id"], json!("ENSG00000141959"));
        assert_eq!(row["cdna_start"], json!(428));
        assert_eq!(row["cdna_end"], json!(428));
        assert_eq!(row["cds_start"], json!(366));
        assert!(row.get("cds_end").is_none(), "{row}");
        // `--protein` and `--symbol` are off, so no core field is renamed to
        // these two keys.
        assert_eq!(row["protein_id"], json!("column"));
        assert_eq!(row["gene_symbol"], json!("column"));
    }

    /// A plugin key named as one of the thirteen columns prints under its own
    /// lowercased name, the columns the object leaves out included
    /// (OutputFactory/JSON.pm:321-343 keeps every key of the consequence hash;
    /// `%SKIP_KEYS` applies to the record's hash alone), and never twice.
    #[test]
    fn a_plugin_key_named_as_a_fixed_column_prints_under_its_name() {
        let formatter = JsonOutputFormatter::new(FieldOptions::default(), Vec::new(), None);
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        tc.plugin_data
            .insert("Location".to_string(), "u".to_string());
        tc.plugin_data
            .insert("Codons".to_string(), "GCC/GTC".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["location"], json!("u"));
        assert_eq!(row["codons"], json!("GCC/GTC"));
        assert_eq!(
            row.as_object()
                .unwrap()
                .keys()
                .filter(|k| *k == "codons")
                .count(),
            1
        );
    }

    /// A column the object leaves out joins it on the rows the plugin wrote the
    /// key alone, never with the core value in the plugin's place: a declared
    /// `Uploaded_variation` or `Location` the plugin left unwritten prints
    /// nothing, as does a `Location` written for the whole variant by a plugin
    /// that runs on regulatory rows alone (OutputFactory.pm:2232), while a
    /// `Uploaded_variation` written on the row prints the plugin's value.
    #[test]
    fn an_omitted_column_a_plugin_left_unwritten_stays_out_of_the_object() {
        let declared = ["Uploaded_variation", "Location"]
            .map(str::to_string)
            .to_vec();
        let formatter = JsonOutputFormatter::new(FieldOptions::default(), declared.clone(), None);
        let v = make_annotated_variant();
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert!(row.get("uploaded_variation").is_none(), "{row}");
        assert!(row.get("location").is_none(), "{row}");

        let regulatory_only = FieldOptions {
            plugin_row_kinds: Arc::new(std::collections::HashMap::from([(
                "Location".to_string(),
                vec![vep_core::consequence::RowKind::RegulatoryFeature],
            )])),
            ..Default::default()
        };
        let formatter = JsonOutputFormatter::new(regulatory_only, declared, None);
        let mut v = make_annotated_variant();
        v.plugin_data
            .insert("Location".to_string(), "u".to_string());
        v.transcript_consequences[0]
            .plugin_data
            .insert("Uploaded_variation".to_string(), "named".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert!(row.get("location").is_none(), "{row}");
        assert_eq!(row["uploaded_variation"], json!("named"));
    }

    /// Ensembl VEP 116.3 writes the AncestralAllele plugin's `AA` as the string
    /// the FASTA index returns for the record's span, empty for a span past the
    /// record's end (AncestralAllele.pm:136-137, :151), and its JSON writer
    /// drops an undefined value alone (OutputFactory/JSON.pm:327), so every
    /// consequence object keeps `"aa": ""`; a record the plugin left without
    /// the key has no key.
    #[test]
    fn an_ancestral_allele_past_the_fasta_record_end_keeps_its_empty_key() {
        let formatter =
            JsonOutputFormatter::new(FieldOptions::default(), vec!["AA".to_string()], None);
        let mut v = make_annotated_variant();
        v.plugin_data.insert("AA".to_string(), String::new());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["aa"], json!(""));
        v.plugin_data.insert("AA".to_string(), "..".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(row["aa"], json!(".."));
        let bare = formatter.format_variant_value(&make_annotated_variant())
            ["transcript_consequences"][0]
            .clone();
        assert!(bare.get("aa").is_none(), "{bare}");
    }

    /// `OutputFactory/JSON.pm:249-258`: `custom_annotations` is a record-level
    /// object keyed by short name, the statistics object unshifted before the
    /// records, each record with `name`, its `allele` when `exact` matched one,
    /// `fields` and `score`, every value numberified (:288) and a `FILTER` list an
    /// array; the record's alleles share one object, and the custom keys appear
    /// in no consequence object, since the factory removes them from the variant
    /// before it builds the entries (:251). A source without a record is listed
    /// for its statistics once another source lists a record (:256), and a
    /// record with no source record writes no object at all.
    #[test]
    fn custom_annotations_is_a_record_level_object_with_stats_first() {
        use crate::output::fields::fixtures::{
            clinvar_annotations, clinvar_source, conservation_annotations, conservation_source,
        };
        let options = FieldOptions {
            custom_sources: vec![clinvar_source(), conservation_source()],
            ..Default::default()
        };
        let mut g = make_annotated_variant();
        g.custom_annotations = vec![clinvar_annotations(b"G"), conservation_annotations()];
        let mut t = make_annotated_variant();
        t.alt_alleles = vec![b"T".to_vec()];
        t.allele_index = 1;
        t.custom_annotations = g.custom_annotations.clone();
        let formatter = JsonOutputFormatter::new(options, Vec::new(), None);
        let value = formatter.format_record_value(&[&g, &t]);
        assert_eq!(
            value["custom_annotations"],
            json!({
                "ClinVar": [
                    {"name": "rs1", "allele": "G", "fields": {"CLNSIG": "Pathogenic", "FILTER": ["q10", "s50"]}},
                    {"name": "21:100-100", "fields": {"CLNSIG": "Benign", "CLNREVSTAT": "criteria_provided"}}
                ],
                "cons": [
                    {"min": 0, "count": 1},
                    {"name": 0, "fields": {"PC": 100}, "score": 0}
                ]
            })
        );
        for entry in value["transcript_consequences"].as_array().unwrap() {
            assert!(entry.get("clinvar").is_none(), "{entry}");
            assert!(entry.get("clinvar_clnsig").is_none(), "{entry}");
            assert!(entry.get("cons_min").is_none(), "{entry}");
        }
        // Statistics of a source without a record ride along; an undefined `PC`
        // and an uncomputed statistics set are `null`.
        let mut sparse = make_annotated_variant();
        let mut cons = conservation_annotations();
        cons.records.clear();
        cons.stats = Some(CustomStats::Undefined);
        let mut clinvar = clinvar_annotations(b"G");
        clinvar.records.truncate(1);
        clinvar.records[0]
            .fields
            .insert("PC".to_string(), CustomFieldValue::Text(String::new()));
        sparse.custom_annotations = vec![cons, clinvar];
        assert_eq!(
            formatter.format_variant_value(&sparse)["custom_annotations"],
            json!({
                "cons": [null],
                "ClinVar": [
                    {"name": "rs1", "allele": "G", "fields": {"CLNSIG": "Pathogenic", "FILTER": ["q10", "s50"], "PC": null}}
                ]
            })
        );
        let mut none = make_annotated_variant();
        let mut stats_only = conservation_annotations();
        stats_only.records.clear();
        none.custom_annotations = vec![stats_only];
        assert!(formatter
            .format_variant_value(&none)
            .get("custom_annotations")
            .is_none());
    }

    /// Ensembl VEP 116.3 pushes each source's records onto its short name's
    /// list (`AnnotationSource/File.pm:437`) and the JSON factory writes that
    /// list whole (`OutputFactory/JSON.pm:249`), so two sources under one
    /// short name that match the same record list it twice, while the
    /// alleles of one input record still share one object.
    #[test]
    fn two_sources_under_one_short_name_each_list_an_equal_record() {
        use crate::output::fields::fixtures::{
            clinvar_annotations, clinvar_source, conservation_annotations,
        };
        let mut second_source = clinvar_source();
        second_source.file = "/data/other.vcf.gz".to_string();
        let options = FieldOptions {
            custom_sources: vec![clinvar_source(), second_source],
            ..Default::default()
        };
        let mut g = make_annotated_variant();
        g.custom_annotations = vec![clinvar_annotations(b"G"), clinvar_annotations(b"G")];
        let mut t = make_annotated_variant();
        t.alt_alleles = vec![b"T".to_vec()];
        t.allele_index = 1;
        t.custom_annotations = g.custom_annotations.clone();
        // A bigWig source annotates each allele by its own span, so one allele
        // may carry an entry the other lacks; it is listed once.
        t.custom_annotations.push(conservation_annotations());
        let formatter = JsonOutputFormatter::new(options, Vec::new(), None);
        let value = formatter.format_record_value(&[&g, &t]);
        let records = value["custom_annotations"]["ClinVar"].as_array().unwrap();
        assert_eq!(records.len(), 4, "{value}");
        assert_eq!(records[0], records[2]);
        assert_eq!(records[1], records[3]);
        assert_eq!(records[0]["name"], "rs1");
        assert_eq!(records[1]["name"], "21:100-100");
        assert_eq!(
            value["custom_annotations"]["cons"],
            json!([{"min": 0, "count": 1}, {"name": 0, "fields": {"PC": 100}, "score": 0}])
        );
    }

    /// A field the source leaves undefined, a GFF window's `feature_id` when
    /// the record has no `ID` attribute (`File/GFF.pm:199`), is `null`: the
    /// JSON factory writes the record's hash whole (`OutputFactory/JSON.pm:249`)
    /// where the text formats join the undefined entry as the empty string.
    #[test]
    fn an_undefined_custom_field_is_null() {
        use crate::output::fields::fixtures::{clinvar_annotations, clinvar_source};
        let options = FieldOptions {
            custom_sources: vec![clinvar_source()],
            ..Default::default()
        };
        let mut v = make_annotated_variant();
        let mut clinvar = clinvar_annotations(b"G");
        clinvar.records[1]
            .fields
            .insert("CLNREVSTAT".to_string(), CustomFieldValue::Undefined);
        v.custom_annotations = vec![clinvar];
        let formatter = JsonOutputFormatter::new(options, Vec::new(), None);
        let value = formatter.format_variant_value(&v);
        assert_eq!(
            value["custom_annotations"]["ClinVar"][1]["fields"],
            json!({"CLNSIG": "Benign", "CLNREVSTAT": null})
        );
    }

    /// `nearest` is a top-level array (`JSON.pm:262` lower-cases the record's
    /// keys at the top level), absent from every consequence object and from a run
    /// without the flag; the empty list is written.
    #[test]
    fn nearest_is_a_top_level_array_and_the_empty_list_is_written() {
        let options = FieldOptions {
            nearest: true,
            ..Default::default()
        };
        let mut v = make_annotated_variant();
        assert!(fmt(options.clone())
            .format_variant_value(&v)
            .get("nearest")
            .is_none());
        v.nearest = Some(vec!["ENST1".into(), "ENST2".into()]);
        let parsed = fmt(options.clone()).format_variant_value(&v);
        assert_eq!(parsed["nearest"], json!(["ENST1", "ENST2"]));
        let tc = &parsed["transcript_consequences"][0];
        assert!(tc.get("nearest").is_none(), "{tc}");
        v.nearest = Some(Vec::new());
        let parsed = fmt(options).format_variant_value(&v);
        assert_eq!(parsed["nearest"], json!([]));
    }

    /// Under `--spdi --ga4gh_vrs` a consequence entry carries `spdi` as a string
    /// and `ga4gh_vrs` as the object itself, both interval bounds numbers, as
    /// Ensembl VEP 116.3 writes with a database for `21 10538885 A C`
    /// (OutputFactory.pm:1409, JSON.pm:288); an allele with no object has no key.
    #[test]
    fn ga4gh_vrs_is_a_nested_object_with_numeric_bounds() {
        let mut v = make_annotated_variant();
        v.spdi = Some("21:25000099:A:G".to_string());
        let options = FieldOptions {
            spdi: true,
            ga4gh_vrs: true,
            ..Default::default()
        };
        let tc = &fmt(options.clone()).format_variant_value(&v)["transcript_consequences"][0];
        assert_eq!(tc["spdi"], "21:25000099:A:G");
        assert!(tc.get("ga4gh_vrs").is_none());
        v.ga4gh_vrs = Some(vep_core::spdi::VrsAllele {
            sequence_id: "NC_000021.9".to_string(),
            start: 25000099,
            end: 25000100,
            sequence: "G".to_string(),
        });
        let tc = &fmt(options).format_variant_value(&v)["transcript_consequences"][0];
        assert_eq!(
            tc["ga4gh_vrs"],
            json!({
                "location": {
                    "type": "SequenceLocation",
                    "interval": {"type": "SimpleInterval", "start": 25000099, "end": 25000100},
                    "sequence_id": "refseq:NC_000021.9"
                },
                "state": {"type": "SequenceState", "sequence": "G"},
                "type": "Allele"
            })
        );
        let plain =
            &fmt(FieldOptions::default()).format_variant_value(&v)["transcript_consequences"][0];
        assert!(plain.get("spdi").is_none() && plain.get("ga4gh_vrs").is_none());
    }

    /// An individual's copy writes `ind` and `zyg` at the top level
    /// (OutputFactory/JSON.pm:261-262 lower-cases the VF-level keys), its own
    /// `allele_string` and `id`; a reference-homozygous copy without
    /// `--process_ref_homs` is an object with no consequence and
    /// `most_severe_consequence` `?`.
    #[test]
    fn individual_copy_writes_ind_and_zyg_at_the_top_level() {
        let mut v = make_annotated_variant();
        v.individual = Some("S2".into());
        v.genotype = vec!["A".into(), "G".into()];
        let options = FieldOptions {
            individual: true,
            ..Default::default()
        };
        let parsed = fmt(options.clone()).format_variant_value(&v);
        assert_eq!(parsed["ind"], "S2");
        assert_eq!(parsed["zyg"], "HET");
        assert_eq!(parsed["id"], "rs123");
        assert!(parsed["transcript_consequences"][0].get("ind").is_none());
        let mut hom_ref =
            InputVariant::without_allele("21".into(), 25000100, "A", std::iter::empty());
        hom_ref.id = Some("rs123".into());
        hom_ref.individual = Some("S3".into());
        hom_ref.genotype = vec!["A".into(), "A".into()];
        hom_ref.hom_ref = true;
        let parsed = fmt(options).format_variant_value(&hom_ref);
        assert_eq!(parsed["ind"], "S3");
        assert_eq!(parsed["zyg"], "HOMREF");
        assert_eq!(parsed["allele_string"], "A");
        assert_eq!(parsed["most_severe_consequence"], "?");
        assert!(parsed.get("transcript_consequences").is_none());
        let plain = fmt(FieldOptions::default()).format_variant_value(&make_annotated_variant());
        assert!(plain.get("ind").is_none() && plain.get("zyg").is_none());
    }

    /// Under `--individual_zyg` the record's `zyg` is the array of `sample:ZYG`
    /// strings (`OutputFactory.pm:981-987` stores a list, which JSON.pm:261-262
    /// writes as is), at the top level and on no consequence entry, with no `ind`;
    /// a record no named sample has a genotype for carries no `zyg`, as Ensembl
    /// VEP sets no key; a record with no allele keeps its object with the list.
    #[test]
    fn individual_zyg_writes_zyg_as_an_array_at_the_top_level() {
        use std::sync::Arc;
        use vep_core::variant::{SampleZygosity, Zygosity, ZygosityList};
        let entry = |sample: &str, zygosity: Zygosity| SampleZygosity {
            sample: sample.into(),
            zygosity,
        };
        let list = |entries: Vec<SampleZygosity>| Some(Arc::new(ZygosityList::new(entries)));
        let mut v = make_annotated_variant();
        v.sample_zygosities = list(vec![
            entry("S1", Zygosity::Het),
            entry("S3", Zygosity::HomRef),
        ]);
        let options = FieldOptions {
            individual_zyg: true,
            ..Default::default()
        };
        let parsed = fmt(options.clone()).format_variant_value(&v);
        assert_eq!(parsed["zyg"], json!(["S1:HET", "S3:HOMREF"]));
        assert!(parsed.get("ind").is_none());
        assert!(parsed["transcript_consequences"][0].get("zyg").is_none());
        v.sample_zygosities = list(Vec::new());
        let parsed = fmt(options.clone()).format_variant_value(&v);
        assert!(parsed.get("zyg").is_none());
        let mut no_allele =
            InputVariant::without_allele("21".into(), 25000100, "A", std::iter::empty());
        no_allele.sample_zygosities = list(vec![entry("S1", Zygosity::HomRef)]);
        let parsed = fmt(options).format_variant_value(&no_allele);
        assert_eq!(parsed["zyg"], json!(["S1:HOMREF"]));
        assert_eq!(parsed["most_severe_consequence"], "?");
    }

    /// A serialised `Zygosity` is its printed label and reads back, so a plugin
    /// handed the serialised record sees the `ZYG` value the output prints.
    #[test]
    fn zygosity_serialises_as_its_printed_label() {
        use vep_core::variant::Zygosity;
        for zygosity in [
            Zygosity::Het,
            Zygosity::Hom,
            Zygosity::HetRef,
            Zygosity::HomRef,
        ] {
            let serialised = serde_json::to_string(&zygosity).unwrap();
            assert_eq!(serialised, format!("\"{}\"", zygosity.label()));
            assert_eq!(
                serde_json::from_str::<Zygosity>(&serialised).unwrap(),
                zygosity
            );
        }
        assert_eq!(
            serde_json::to_string(&Zygosity::HomRef).unwrap(),
            "\"HOMREF\""
        );
    }

    /// An unannotated record keeps VEP's input-only shape plus the identity key.
    #[test]
    fn skipped_record_is_its_input_line_and_the_version() {
        let mut v = make_annotated_variant();
        v.vep_skip = true;
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        let keys: Vec<&String> = parsed.as_object().unwrap().keys().collect();
        assert_eq!(keys, ["input", "vep_rs_version"]);
        assert_eq!(parsed["vep_rs_version"], fields::VEP_RS_VERSION);
    }

    /// A record with no allele has no first allele to read the identity from:
    /// its value is the object carrying the version alone, whatever assembly
    /// the formatter names.
    #[test]
    fn a_record_with_no_allele_is_the_version_object() {
        let value = fmt(FieldOptions::default()).format_record_value(&[]);
        assert_eq!(value, json!({ "vep_rs_version": fields::VEP_RS_VERSION }));
    }

    /// The entries that write lines write none for a record with no allele.
    #[test]
    fn a_record_with_no_allele_writes_no_line() {
        let formatter = fmt(FieldOptions::default());
        assert_eq!(formatter.format_record(&[]).unwrap(), Vec::<String>::new());
        assert_eq!(
            formatter.format_record_bytes(&[]).unwrap(),
            Vec::<Vec<u8>>::new()
        );
    }

    #[test]
    fn transcript_consequence_uses_vep_keys_and_flag_gating() {
        let v = make_annotated_variant();
        let plain = fmt(FieldOptions::default()).format_variant_value(&v);
        let tc = &plain["transcript_consequences"][0];
        assert_eq!(tc["variant_allele"], "G");
        assert_eq!(tc["consequence_terms"][0], "synonymous_variant");
        assert_eq!(tc["impact"], "LOW");
        assert_eq!(tc["gene_id"], "ENSG00000141959");
        assert_eq!(tc["transcript_id"], "ENST00000352957");
        assert_eq!(tc["strand"], -1);
        assert_eq!(tc["cdna_start"], 428);
        assert_eq!(tc["cdna_end"], 428);
        assert_eq!(tc["cds_start"], 366);
        assert_eq!(tc["protein_start"], 122);
        assert_eq!(tc["amino_acids"], "L");
        assert_eq!(tc["codons"], "ctG/ctA");
        assert_eq!(tc["flags"][0], "cds_end_NF");
        assert!(tc.get("feature_type").is_none());
        assert!(tc.get("gene_symbol").is_none(), "SYMBOL needs --symbol");
        assert!(tc.get("biotype").is_none(), "BIOTYPE needs --biotype");
        assert!(tc.get("canonical").is_none(), "CANONICAL needs --canonical");

        let flagged = fmt(FieldOptions {
            symbol: true,
            biotype: true,
            canonical: true,
            ..Default::default()
        })
        .format_variant_value(&v);
        let tc = &flagged["transcript_consequences"][0];
        assert_eq!(tc["gene_symbol"], "MRPL39");
        assert_eq!(tc["gene_symbol_source"], "HGNC");
        assert_eq!(tc["hgnc_id"], "HGNC:14027");
        assert_eq!(tc["biotype"], "protein_coding");
        assert_eq!(tc["canonical"], 1);
    }

    #[test]
    fn intergenic_record_has_an_intergenic_consequences_array() {
        let mut v = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        assert_eq!(parsed["most_severe_consequence"], "intergenic_variant");
        assert!(parsed.get("transcript_consequences").is_none());
        let ic = &parsed["intergenic_consequences"][0];
        assert_eq!(ic["consequence_terms"][0], "intergenic_variant");
        assert_eq!(ic["impact"], "MODIFIER");
        assert_eq!(ic["variant_allele"], "G");
    }

    /// `--shift_genomic 1` puts a `location` on every transcript entry, moved by
    /// the row's shift along the transcript's strand; an intergenic entry gets
    /// one only when `--shift_3prime 1` analysed its allele as an insertion or
    /// deletion, and nothing gets one without `--shift_genomic 1`.
    #[test]
    fn shift_genomic_writes_a_location_on_the_entries_vep_moves() {
        let shifting = |shift_3prime: bool| FieldOptions {
            shift_3prime,
            shift_genomic: true,
            ..Default::default()
        };
        let mut v = make_annotated_variant();
        v.transcript_consequences[0].shift_length = 2;
        let plain = fmt(FieldOptions::default()).format_variant_value(&v);
        assert!(plain["transcript_consequences"][0]
            .get("location")
            .is_none());
        let moved = fmt(shifting(true)).format_variant_value(&v);
        assert_eq!(
            moved["transcript_consequences"][0]["location"],
            "21:25000098"
        );
        let unmoved = fmt(shifting(false)).format_variant_value(&v);
        assert_eq!(
            unmoved["transcript_consequences"][0]["location"],
            "21:25000100"
        );

        let mut intergenic = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"-".to_vec(),
        );
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let unshifted = fmt(shifting(true)).format_variant_value(&intergenic);
        assert!(unshifted["intergenic_consequences"][0]
            .get("location")
            .is_none());
        intergenic.intergenic_shift = Some(3);
        let shifted = fmt(shifting(true)).format_variant_value(&intergenic);
        assert_eq!(
            shifted["intergenic_consequences"][0]["location"],
            "21:25000103"
        );
        let still = fmt(shifting(false)).format_variant_value(&intergenic);
        assert!(still["intergenic_consequences"][0]
            .get("location")
            .is_none());
    }

    /// A structural variant's transcript entries are built by
    /// `TranscriptStructuralVariationAllele_to_output_hash`
    /// (OutputFactory.pm:2127-2134), which writes no `Location`, so
    /// `--shift_genomic 1` puts no `location` on them whether or not
    /// `--shift_3prime 1` is on.
    #[test]
    fn shift_genomic_writes_no_location_on_a_structural_variants_entries() {
        let mut sv = make_annotated_variant();
        sv.is_structural = true;
        sv.end = 25000400;
        sv.sv_end = Some(25000400);
        sv.alt_alleles = vec![b"<DEL>".to_vec()];
        sv.record_so_term = Some("deletion".to_string());
        sv.transcript_consequences[0].consequences = smallvec![Consequence::FeatureTruncation];
        for shift_3prime in [false, true] {
            let options = FieldOptions {
                shift_3prime,
                shift_genomic: true,
                ..Default::default()
            };
            let value = fmt(options).format_variant_value(&sv);
            assert!(
                value["transcript_consequences"][0]
                    .get("location")
                    .is_none(),
                "shift_3prime {shift_3prime}: {value}"
            );
        }
    }

    /// A `--summary` or `--most_severe` record files one entry under
    /// `intergenic_consequences` holding `consequence_terms` alone
    /// (`add_VariationFeatureOverlapAllele_info` hands `summary_only` an empty
    /// hash, JSON.pm:320 and 348), its most severe consequence is the first term,
    /// and the record keeps its variant-level keys and co-located variants.
    #[test]
    fn summary_entry_is_intergenic_consequences_with_consequence_terms_only() {
        let mut v = make_annotated_variant();
        v.transcript_consequences.clear();
        v.most_severe_consequence = Some(Consequence::SynonymousVariant);
        v.summary_terms = Some(smallvec![
            Consequence::SynonymousVariant,
            Consequence::IntronVariant
        ]);
        v.colocated_variants = vec![{
            let mut colocated = vep_core::variant::ColocatedVariant::default();
            colocated.id = "rs123".into();
            colocated
        }];
        let options = FieldOptions {
            variant_class: true,
            ..Default::default()
        };
        let parsed = fmt(options).format_variant_value(&v);
        assert!(parsed.get("transcript_consequences").is_none());
        let entries = parsed["intergenic_consequences"].as_array().unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(
            entries[0],
            json!({"consequence_terms": ["synonymous_variant", "intron_variant"]})
        );
        assert_eq!(parsed["most_severe_consequence"], "synonymous_variant");
        assert_eq!(parsed["variant_class"], "SNV");
        assert_eq!(parsed["colocated_variants"][0]["id"], "rs123");
    }

    #[test]
    fn multi_allelic_record_is_one_object() {
        let mut a = make_annotated_variant();
        a.id = Some(".".into());
        a.allele_index = 0;
        a.uploaded_allele_string = Some("A/G/T".into());
        a.record_allele_string_multi = Some("A/G/T".into());
        let mut b = a.clone();
        b.alt_alleles = vec![b"T".to_vec()];
        b.allele_string = "A/T".into();
        b.allele_index = 1;
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert_eq!(parsed["allele_string"], "A/G/T");
        // A VCF record whose ID column is `.` keeps the dot as its JSON id.
        assert_eq!(parsed["id"], ".");
        let tcs = parsed["transcript_consequences"].as_array().unwrap();
        assert_eq!(tcs.len(), 2);
        assert_eq!(tcs[0]["variant_allele"], "G");
        assert_eq!(tcs[1]["variant_allele"], "T");
    }

    /// The consequence object of a variant whose three position columns read
    /// `cdna`, `cds` and `protein`.
    fn positions(cdna: &str, cds: &str, protein: &str) -> Value {
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        tc.cdna_position = Some(cdna.into());
        tc.cds_position = Some(cds.into());
        tc.protein_position = Some(protein.into());
        fmt(FieldOptions::default()).format_variant_value(&v)["transcript_consequences"][0].clone()
    }

    /// An undefined end is left out: `445-?`, `172-?` and `58-?` give
    /// `cdna_start` 445, `cds_start` 172 and `protein_start` 58 with no `*_end`
    /// key. The value follows ensembl-vep `OutputFactory/JSON.pm:337-341`
    /// (release 116.3), which keeps a `?` end as the end and then deletes any end
    /// that is not a number, and the row of the same shape Ensembl VEP 116.3
    /// writes for `tests/golden/116/GRCh38` record `2561884`
    /// (`1 15443553 GGAGT G`) against ENST00000375949: `514-?`, `492-?` and
    /// `164-?` give `cdna_start` 514, `cds_start` 492 and `protein_start` 164
    /// alone.
    #[test]
    fn undefined_end_writes_the_start_alone() {
        let tc = positions("445-?", "172-?", "58-?");
        assert_eq!(tc["cdna_start"], 445);
        assert_eq!(tc["cds_start"], 172);
        assert_eq!(tc["protein_start"], 58);
        for key in ["cdna_end", "cds_end", "protein_end"] {
            assert!(
                tc.get(key).is_none(),
                "{key} is written for an undefined end"
            );
        }
    }

    /// An undefined start is left out: `?-13`, `?-13` and `?-5` give `cdna_end`
    /// 13, `cds_end` 13 and `protein_end` 5 with no `*_start` key. The value
    /// follows ensembl-vep `OutputFactory/JSON.pm:337-341` (release 116.3), which
    /// deletes any start that is not a number, and the row of the same shape
    /// Ensembl VEP 116.3 writes for `tests/golden/116/GRCh38` record `4287535`
    /// (`1 10434655 AGAAAATTTTGCCAAAGACCTT A`) against ENST00000309048: `?-280`,
    /// `?-195` and `?-65` give `cdna_end` 280, `cds_end` 195 and `protein_end` 65
    /// alone.
    #[test]
    fn undefined_start_writes_the_end_alone() {
        let tc = positions("?-13", "?-13", "?-5");
        assert_eq!(tc["cdna_end"], 13);
        assert_eq!(tc["cds_end"], 13);
        assert_eq!(tc["protein_end"], 5);
        for key in ["cdna_start", "cds_start", "protein_start"] {
            assert!(
                tc.get(key).is_none(),
                "{key} is written for an undefined start"
            );
        }
    }

    /// `OutputFactory/JSON.pm`: a defined range writes both ends.
    #[test]
    fn defined_range_writes_both_ends() {
        let tc = positions("1161-1162", "900-901", "300-301");
        assert_eq!(
            (&tc["cdna_start"], &tc["cdna_end"]),
            (&json!(1161), &json!(1162))
        );
        assert_eq!(
            (&tc["cds_start"], &tc["cds_end"]),
            (&json!(900), &json!(901))
        );
        assert_eq!(
            (&tc["protein_start"], &tc["protein_end"]),
            (&json!(300), &json!(301))
        );
    }

    /// `OutputFactory/JSON.pm`: a single position is both start and end.
    #[test]
    fn single_position_writes_both_ends() {
        let tc = positions("428", "366", "122");
        assert_eq!(
            (&tc["cdna_start"], &tc["cdna_end"]),
            (&json!(428), &json!(428))
        );
        assert_eq!(
            (&tc["cds_start"], &tc["cds_end"]),
            (&json!(366), &json!(366))
        );
        assert_eq!(
            (&tc["protein_start"], &tc["protein_end"]),
            (&json!(122), &json!(122))
        );
    }

    /// The `--total_length` forms. `N/T` is kept as a string at both ends and
    /// `N-?/T` writes the start alone, as `OutputFactory/JSON.pm` writes them
    /// (`numberify` leaves a string with a slash alone). A range keeps its end:
    /// `N-M/T` writes `N` and `M/T`, and `?-N/T` writes `N/T` under the end key
    /// alone, where Ensembl VEP writes `N` at both ends of the first and neither
    /// key of the second (`JSON.pm:337` takes the start in place of an end that
    /// carries a total).
    #[test]
    fn total_length_forms_keep_the_end_of_a_range() {
        assert_eq!(
            position_bounds("428/1234"),
            (Some(json!("428/1234")), Some(json!("428/1234")))
        );
        assert_eq!(
            position_bounds("428-430/1234"),
            (Some(json!(428)), Some(json!("430/1234")))
        );
        assert_eq!(position_bounds("428-?/1234"), (Some(json!(428)), None));
        assert_eq!(
            position_bounds("?-430/1234"),
            (None, Some(json!("430/1234")))
        );
        assert_eq!(position_bounds("-/1234"), (None, None));
    }

    /// The `--total_length` columns as the JSON writer files them: a single
    /// position with its total is the string `"53/500"` at both ends, Ensembl
    /// VEP's shape (`OutputFactory/JSON.pm`
    /// `add_VariationFeatureOverlapAllele_info`); a range writes its start as a
    /// number and its end with the total, `"6/300"`, where Ensembl VEP writes the
    /// start under both keys; `consequence_terms` is untouched.
    #[test]
    fn total_length_json_writes_a_range_end_with_its_total() {
        let mut v = make_annotated_variant();
        let consequence = &mut v.transcript_consequences[0];
        consequence.cdna_position = Some("53".into());
        consequence.cds_position = Some("4-6".into());
        consequence.protein_position = Some("2".into());
        consequence.cdna_total = Some(500);
        consequence.cds_total = Some(300);
        consequence.protein_total = Some(99);
        let tc = fmt(FieldOptions {
            total_length: true,
            ..Default::default()
        })
        .format_variant_value(&v)["transcript_consequences"][0]
            .clone();
        assert_eq!(
            (&tc["cdna_start"], &tc["cdna_end"]),
            (&json!("53/500"), &json!("53/500"))
        );
        assert_eq!(
            (&tc["cds_start"], &tc["cds_end"]),
            (&json!(4), &json!("6/300"))
        );
        assert_eq!(
            (&tc["protein_start"], &tc["protein_end"]),
            (&json!("2/99"), &json!("2/99"))
        );
        assert_eq!(tc["consequence_terms"], json!(["synonymous_variant"]));
    }

    /// `--terms display` renders `consequence_terms`, `most_severe_consequence`
    /// and the intergenic entry in the vocabulary; `variant_class` stays the SO
    /// class term, as `class_SO_term` is what Ensembl VEP writes there.
    #[test]
    fn terms_vocabulary_renders_json_consequences() {
        let mut v = make_annotated_variant();
        v.transcript_consequences[0].consequences = smallvec![
            Consequence::SpliceRegionVariant,
            Consequence::SynonymousVariant
        ];
        v.most_severe_consequence = Some(Consequence::SpliceRegionVariant);
        let parsed = fmt(FieldOptions {
            variant_class: true,
            terms: Vocabulary::Display,
            ..Default::default()
        })
        .format_variant_value(&v);
        assert_eq!(
            parsed["transcript_consequences"][0]["consequence_terms"],
            json!(["SPLICE_SITE", "SYNONYMOUS_CODING"])
        );
        assert_eq!(parsed["most_severe_consequence"], "SPLICE_SITE");
        assert_eq!(parsed["variant_class"], "SNV");
        let mut intergenic = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let parsed = fmt(FieldOptions {
            terms: Vocabulary::Ncbi,
            ..Default::default()
        })
        .format_variant_value(&intergenic);
        assert_eq!(
            parsed["intergenic_consequences"][0]["consequence_terms"],
            json!(["intergenic_variant"])
        );
        assert_eq!(parsed["most_severe_consequence"], "intergenic_variant");
    }

    /// `MINIMISED` and `AMBIGUITY`, like `VARIANT_CLASS`, are the variant's own
    /// keys (`OutputFactory.pm` `VariationFeature_to_output_hash`), so the JSON
    /// writer puts them on the record object, lower-cased and numberified, and on
    /// no consequence entry (`OutputFactory/JSON.pm` copies the variant hash onto
    /// the record and builds the entries from an empty hash).
    #[test]
    fn variant_level_keys_sit_on_the_record_object() {
        let mut v = make_annotated_variant();
        v.minimised = true;
        let parsed = fmt(FieldOptions {
            minimal: true,
            ambiguity: true,
            dont_skip: true,
            ..Default::default()
        })
        .format_variant_value(&v);
        assert_eq!(parsed["minimised"], json!(1));
        assert_eq!(parsed["ambiguity"], json!("R"));
        assert!(parsed.get("check_ref").is_none());
        let entry = &parsed["transcript_consequences"][0];
        for key in ["minimised", "ambiguity", "check_ref", "variant_class"] {
            assert!(entry.get(key).is_none(), "{key} on a consequence entry");
        }
    }

    /// A `colocated_variants` entry as `add_colocated_variant_info_JSON` writes
    /// it (OutputFactory/JSON.pm:407-476): `var_synonyms` nested by source with
    /// its ids split on `,` and numberified where they look like numbers,
    /// `clin_sig` and `pubmed` as lists, `clin_sig_allele` and
    /// `clin_sig_ref_allele` as written, `clinical_impact` as the raw string
    /// under `clinvar_somatic_classification`, `chr` as `seq_region_name`, the
    /// zero flags left out, whatever the flags; the consequence objects carry none
    /// of the text formats' co-located keys (:493-515), and the record carries no
    /// `existing_variation`. Ensembl VEP 116.3 writes
    /// `{"var_synonyms": {"ClinVar": ["RCV000886041", "VCV000713837"]}, "clin_sig":
    /// ["likely_benign"], "clin_sig_allele": "T:likely_benign", "clin_sig_ref_allele":
    /// "G", "phenotype_or_disease": 1, "id": "rs374325784", "seq_region_name": "21",
    /// "start": 10541121, "end": 10541121, "strand": 1, "allele_string": "G/A/C/T"}`
    /// for `21 10541121 G A` of the `GRCh38-colocated-alt` corpus.
    #[test]
    fn var_synonyms_nested_by_source() {
        let mut v = make_annotated_variant();
        v.colocated_variants = vec![{
            let mut colocated = vep_core::variant::ColocatedVariant::default();
            colocated.chr = Some("21".into());
            colocated.id = "rs374325784".into();
            colocated.start = 10_541_121;
            colocated.end = 10_541_121;
            colocated.allele_string = Some("G/A/C/T".into());
            colocated.strand = 1;
            colocated.clin_sig = Some("likely_benign".into());
            colocated.phenotype = true;
            colocated.pubmed = Some("26076356,11524731".into());
            colocated.clinvar = Some(Box::new({
                let mut detail = vep_core::variation::ClinVarDetail::default();
                detail.clin_sig_allele = Some("T:likely_benign".into());
                detail.clin_sig_ref_allele = Some("G".into());
                detail.clinical_impact = Some("Neoplasm(Tier I - Strong)".into());
                detail
            }));
            colocated.var_synonyms = Some(
                "ClinVar::RCV000886041,VCV000713837--OMIM::605717.0001--UniProt::VAR_089699".into(),
            );
            colocated
        }];
        v.existing_variation = vec!["rs374325784".into()];
        let options = FieldOptions {
            check_existing: true,
            pubmed: true,
            ..Default::default()
        };
        let parsed = fmt(options.clone()).format_variant_value(&v);
        assert_eq!(
            parsed["colocated_variants"][0],
            json!({
                "seq_region_name": "21",
                "id": "rs374325784",
                "start": 10541121,
                "end": 10541121,
                "strand": 1,
                "allele_string": "G/A/C/T",
                "phenotype_or_disease": 1,
                "clin_sig": ["likely_benign"],
                "pubmed": [26076356, 11524731],
                "clin_sig_allele": "T:likely_benign",
                "clin_sig_ref_allele": "G",
                "clinvar_somatic_classification": "Neoplasm(Tier I - Strong)",
                "var_synonyms": {
                    "ClinVar": ["RCV000886041", "VCV000713837"],
                    "OMIM": [605717.0001],
                    "UniProt": ["VAR_089699"]
                }
            })
        );
        assert!(parsed.get("existing_variation").is_none());
        let tc = &parsed["transcript_consequences"][0];
        for key in [
            "clin_sig",
            "pheno",
            "somatic",
            "pubmed",
            "existing_variation",
        ] {
            assert!(tc.get(key).is_none(), "{key} in a consequence object");
        }
        // A stored list is split as `split(',', ...)` splits it (:469): no empty
        // last element, `[]` for separators alone, no key for a false value.
        let lists = |pubmed: &str, clin_sig: &str| {
            let mut stored = v.clone();
            stored.colocated_variants[0].pubmed = Some(pubmed.into());
            stored.colocated_variants[0].clin_sig = Some(clin_sig.into());
            let entry =
                fmt(options.clone()).format_variant_value(&stored)["colocated_variants"][0].clone();
            (entry.get("pubmed").cloned(), entry.get("clin_sig").cloned())
        };
        assert_eq!(
            lists("99000001,99000002,", "benign,likely_benign,"),
            (
                Some(json!([99000001, 99000002])),
                Some(json!(["benign", "likely_benign"]))
            )
        );
        assert_eq!(
            lists(",99000001", "benign,,pathogenic"),
            (
                Some(json!(["", 99000001])),
                Some(json!(["benign", "", "pathogenic"]))
            )
        );
        assert_eq!(lists(",,", ","), (Some(json!([])), Some(json!([]))));
        assert_eq!(lists("0", "0"), (None, None));
        assert_eq!(
            lists("0,0", ",,0"),
            (Some(json!([0, 0])), Some(json!(["", "", 0])))
        );
        // A germline record without ClinVar fields or synonyms carries the
        // coordinates, alleles and name alone (:433-435 removes the rest).
        v.colocated_variants = vec![{
            let mut colocated = vep_core::variant::ColocatedVariant::default();
            colocated.chr = Some("21".into());
            colocated.id = "rs1300916113".into();
            colocated.start = 5_030_602;
            colocated.end = 5_030_602;
            colocated.allele_string = Some("A/C".into());
            colocated.strand = 1;
            colocated
        }];
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        assert_eq!(
            parsed["colocated_variants"][0],
            json!({
                "seq_region_name": "21",
                "id": "rs1300916113",
                "start": 5030602,
                "end": 5030602,
                "strand": 1,
                "allele_string": "A/C"
            })
        );
        assert_eq!(
            var_synonyms_by_source("A::1,2--B--A::3--"),
            json!({"A": [3], "B": []}),
            "a repeated source keeps its last group, a group without `::` holds no id, a trailing empty group is no group"
        );
    }

    /// A cache column outside the fixed set and the frequency populations is
    /// echoed under its own name, numberified: `add_colocated_variant_info_JSON`
    /// copies the record's hash (OutputFactory/JSON.pm:415), prunes an empty
    /// value and a `0` under a key without `af` in it (:433-435) and removes the
    /// `@FREQ_FIELDS` populations it nested under `frequencies` (:439-441).
    /// Ensembl VEP 116.3 writes `minor_allele` and `minor_allele_freq` in every
    /// object on a cache declaring the two pre-108 columns; a release-116 cache
    /// declares no such column.
    #[test]
    fn colocated_variant_echoes_the_other_cache_columns() {
        let mut v = make_annotated_variant();
        v.colocated_variants = vec![{
            let mut colocated = vep_core::variant::ColocatedVariant::default();
            colocated.chr = Some("21".into());
            colocated.id = "rs1".into();
            colocated.start = 1_004;
            colocated.end = 1_004;
            colocated.allele_string = Some("C/A".into());
            colocated.strand = 1;
            colocated.frequency_strings = [
                ("minor_allele", "A"),
                ("minor_allele_freq", "0.12"),
                ("AF", "A:0.5"),
                ("other_af", "0"),
                ("some_count", "0"),
                ("some_text", "x"),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
            colocated
        }];
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        let entry = &parsed["colocated_variants"][0];
        assert_eq!(entry["minor_allele"], "A");
        assert_eq!(entry["minor_allele_freq"], 0.12);
        assert_eq!(entry["other_af"], 0, "a `0` under a key with `af` stays");
        assert_eq!(entry["some_text"], "x");
        assert!(entry.get("some_count").is_none(), "a `0` is pruned");
        assert!(
            entry.get("AF").is_none(),
            "a population is written under frequencies, never echoed"
        );
    }

    /// `frequencies` nests, under each alternate allele of the record, the
    /// values of the record allele matched to it keyed by the lower-cased
    /// population (OutputFactory/JSON.pm:419-431), whatever the flags
    /// (:160-163): Ensembl VEP 116.3 writes
    /// `{"A": {"gnomadg": 0.0002047, "gnomadg_afr": 0, ...}, "G": {"af": 1,
    /// "afr": 1, ...}}` for `rs1555870136` on `21 7759128 T A,G` of the
    /// `GRCh38-colocated` corpus with no frequency flag. `AF` prints with four
    /// decimals and numberified (`0.006789` is `0.0068`,
    /// OutputFactory.pm:1219; `rs1555830109` on `21 5131441 T C`); an allele
    /// the variant was not matched to has no entry, and a record without a
    /// frequency carries no `frequencies` key.
    #[test]
    fn colocated_frequencies_nest_per_allele() {
        let mut a = make_annotated_variant();
        a.allele_index = 0;
        a.record_allele_string_multi = Some("A/G/T".into());
        let matched = |input: &str, index: usize| vep_core::variant::MatchedAllele {
            input_allele: input.to_string(),
            input_index: index,
            record_allele: input.to_string(),
            record_index: index,
        };
        let mut rs1 = vep_core::variant::ColocatedVariant::default();
        rs1.chr = Some("21".into());
        rs1.id = "rs1".into();
        rs1.start = 25000100;
        rs1.end = 25000100;
        rs1.allele_string = Some("A/G/T".into());
        rs1.strand = 1;
        rs1.matched_alleles = vec![matched("G", 0), matched("T", 1)];
        rs1.frequency_strings = [
            ("AF", "G:1"),
            ("gnomADg_AFR", "T:0,G:1"),
            ("gnomADg", "T:0.0002047,G:0.9994"),
        ]
        .into_iter()
        .map(|(column, value)| (column.to_string(), value.to_string()))
        .collect();
        a.colocated_variants = vec![rs1.clone()];
        let mut b = a.clone();
        b.alt_alleles = vec![b"T".to_vec()];
        b.allele_string = "A/T".into();
        b.allele_index = 1;
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert_eq!(
            parsed["colocated_variants"][0]["frequencies"],
            json!({
                "G": {"af": 1, "gnomadg": 0.9994, "gnomadg_afr": 1},
                "T": {"gnomadg": 0.0002047, "gnomadg_afr": 0}
            })
        );

        let mut rounded = rs1.clone();
        rounded.matched_alleles = vec![matched("G", 0)];
        rounded.frequency_strings = [("AF".to_string(), "G:0.006789".to_string())].into();
        a.colocated_variants = vec![rounded.clone()];
        b.colocated_variants = vec![rounded];
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert_eq!(
            parsed["colocated_variants"][0]["frequencies"],
            json!({"G": {"af": 0.0068}})
        );

        let mut none = rs1;
        none.frequency_strings.clear();
        a.colocated_variants = vec![none.clone()];
        b.colocated_variants = vec![none];
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert!(parsed["colocated_variants"][0].get("frequencies").is_none());
    }

    /// The alleles the frequencies are read for are those of the record's own
    /// allele string (OutputFactory/JSON.pm:276-284 splits `allele_string`), so
    /// under `--minimal` a split record's trimmed allele is looked up by its
    /// written spelling and finds no entry: Ensembl VEP 116.3 writes
    /// rs1370167485 on `21 5030602 A AA,G` of the `GRCh38-colocated` corpus with
    /// no `frequencies` key under `--check_existing --minimal`, and with
    /// `{"A": {...}}` on the bi-allelic `21 5030602 A AA`, whose string is `-/A`.
    #[test]
    fn colocated_frequencies_are_keyed_by_the_records_own_alleles() {
        let matched = vep_core::variant::MatchedAllele {
            input_allele: "A".into(),
            input_index: 0,
            record_allele: "AA".into(),
            record_index: 0,
        };
        let mut rs = vep_core::variant::ColocatedVariant::default();
        rs.chr = Some("21".into());
        rs.id = "rs1370167485".into();
        rs.start = 5030602;
        rs.end = 5030602;
        rs.allele_string = Some("A/AA".into());
        rs.strand = 1;
        rs.matched_alleles = vec![matched];
        rs.frequency_strings = [("gnomADg".to_string(), "AA:0.0001735".to_string())].into();
        let mut split =
            InputVariant::new("21".into(), 5030603, 5030602, b"-".to_vec(), b"A".to_vec());
        split.record_allele_string_multi = Some("A/AA/G".into());
        split.original_start = Some(5030602);
        split.original_end = Some(5030602);
        split.minimised = true;
        split.colocated_variants = vec![rs.clone()];
        let mut g = split.clone();
        g.allele_index = 1;
        g.alt_alleles = vec![b"G".to_vec()];
        g.allele_string = "A/G".into();
        g.start = 5030602;
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&split, &g]);
        assert_eq!(parsed["allele_string"], "A/AA/G");
        assert!(parsed["colocated_variants"][0].get("frequencies").is_none());

        let mut bi_allelic =
            InputVariant::new("21".into(), 5030603, 5030602, b"-".to_vec(), b"A".to_vec());
        bi_allelic.minimised = true;
        bi_allelic.colocated_variants = vec![rs];
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&bi_allelic]);
        assert_eq!(parsed["allele_string"], "-/A");
        assert_eq!(
            parsed["colocated_variants"][0]["frequencies"],
            json!({"A": {"gnomadg": 0.0001735}})
        );
    }

    /// The consequence object of the annotated variant under the attribute
    /// flags, with the given cross-reference values on its transcript.
    fn attribute_object(edit: impl Fn(&mut TranscriptConsequence)) -> Value {
        let mut v = make_annotated_variant();
        edit(&mut v.transcript_consequences[0]);
        let options = FieldOptions {
            uniprot: true,
            xref_refseq: true,
            mirna: true,
            mane: true,
            mane_select: true,
            mane_plus_clinical: true,
            ..Default::default()
        };
        fmt(options).format_variant_value(&v)["transcript_consequences"][0].clone()
    }

    /// `SWISSPROT`, `TREMBL`, `UNIPARC` and `UNIPROT_ISOFORM` are arrays in
    /// Ensembl VEP's JSON (`[split(',', $id)]`, OutputFactory.pm:1581 at
    /// release/116.3, encoded as is by `JSON.pm`), one element per
    /// comma-separated accession, and a `-` value sets no key: the GRCh37 APP
    /// transcript's `TREMBL=L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN` is a
    /// three-element array.
    #[test]
    fn uniprot_columns_are_arrays() {
        let tc = attribute_object(|tc| {
            edit_transcript(tc, |t| {
                t.swissprot = Some("A4_HUMAN".to_string());
                t.trembl = Some("L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN".to_string());
                t.uniparc = Some("UPI000002DB1C".to_string());
                t.uniprot_isoform = Some("-".to_string());
            });
        });
        assert_eq!(tc["swissprot"], json!(["A4_HUMAN"]));
        assert_eq!(
            tc["trembl"],
            json!(["L7XE61_HUMAN", "L7XCZ9_HUMAN", "B7Z313_HUMAN"])
        );
        assert_eq!(tc["uniparc"], json!(["UPI000002DB1C"]));
        assert!(tc.get("uniprot_isoform").is_none(), "{tc}");
        let empty_element =
            attribute_object(|tc| edit_transcript(tc, |t| t.trembl = Some("A,,B".to_string())));
        assert_eq!(empty_element["trembl"], json!(["A", "", "B"]));
    }

    /// `RefSeq` is renamed `refseq_transcript_ids` (JSON.pm:106) and is the
    /// array of the comma-separated ids (OutputFactory.pm:1549):
    /// `RefSeq=NM_198996.4,NM_001302998.2` on the LIPI transcript
    /// ENST00000681601 is a two-element array.
    #[test]
    fn refseq_transcript_ids_is_array() {
        let tc = attribute_object(|tc| {
            edit_transcript(tc, |t| {
                t.refseq = Some("NM_198996.4,NM_001302998.2".to_string())
            });
        });
        assert_eq!(
            tc["refseq_transcript_ids"],
            json!(["NM_198996.4", "NM_001302998.2"])
        );
        assert!(tc.get("refseq").is_none());
        let dash = attribute_object(|tc| edit_transcript(tc, |t| t.refseq = Some("-".to_string())));
        assert!(dash.get("refseq_transcript_ids").is_none());
    }

    /// The RefSeq and merged cache keys are the lower-cased column names
    /// (`$vfoa_hash->{lc($key)}`, JSON.pm:323-345): `refseq_match` keeps the
    /// list `OutputFactory.pm:1557` builds as an array, `source` and `bam_edit`
    /// are the column texts, and `refseq_offset` is a number, as `numberify`
    /// leaves the integer the hash holds; a transcript without a code, source,
    /// status or offset sets no key, and the keys follow the flags as the text
    /// columns do.
    #[test]
    fn refseq_match_is_an_array_and_the_refseq_keys_are_lower_case() {
        use vep_core::consequence::SourceCache;
        use vep_core::transcript::BamEditStatus;
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        edit_transcript(tc, |t| {
            t.attributes = with_rseq_codes(&["rseq_5p_mismatch", "rseq_cds_mismatch"]).attributes;
            t.source_cache = Some(SourceCache::RefSeq);
        });
        tc.bam_edit_status = Some(BamEditStatus::Ok);
        tc.refseq_offset = Some(-3);
        let merged = FieldOptions {
            merged: true,
            bam_edited: true,
            ..Default::default()
        };
        let object = fmt(merged).format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            object["refseq_match"],
            json!(["rseq_5p_mismatch", "rseq_cds_mismatch"])
        );
        assert_eq!(object["source"], "RefSeq");
        assert_eq!(object["bam_edit"], "OK");
        assert_eq!(object["refseq_offset"], json!(-3));
        for upper in ["REFSEQ_MATCH", "SOURCE", "BAM_EDIT", "REFSEQ_OFFSET"] {
            assert!(object.get(upper).is_none(), "{upper}");
        }

        let refseq = FieldOptions {
            refseq: true,
            ..Default::default()
        };
        let object = fmt(refseq).format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            object["refseq_match"],
            json!(["rseq_5p_mismatch", "rseq_cds_mismatch"])
        );
        assert_eq!(object["refseq_offset"], json!(-3));
        assert!(object.get("source").is_none(), "{object}");
        // `BAM_EDIT` is set without a flag check, so the key survives the group.
        assert_eq!(object["bam_edit"], "OK");

        let plain = make_annotated_variant();
        let object = fmt(FieldOptions {
            merged: true,
            ..Default::default()
        })
        .format_variant_value(&plain)["transcript_consequences"][0]
            .clone();
        for key in ["refseq_match", "source", "bam_edit", "refseq_offset"] {
            assert!(object.get(key).is_none(), "{key}: {object}");
        }
    }

    /// `mirna` is the array of structure terms (OutputFactory.pm:1671), `mane`
    /// the array of set names (:1601, :1609); an undefined element is `null`.
    #[test]
    fn mirna_and_mane_are_arrays() {
        let tc = attribute_object(|tc| {
            tc.mirna = vec!["miRNA_loop".to_string(), "miRNA_stem".to_string()];
            edit_transcript(tc, |t| {
                t.mane_select = Some("NM_021075.4".to_string());
                t.mane_plus_clinical = Some("NM_058174.3".to_string());
            });
        });
        assert_eq!(tc["mirna"], json!(["miRNA_loop", "miRNA_stem"]));
        assert_eq!(tc["mane"], json!(["MANE_Select", "MANE_Plus_Clinical"]));
        assert_eq!(tc["mane_select"], "NM_021075.4");
        assert_eq!(tc["mane_plus_clinical"], "NM_058174.3");
        let past_the_end = attribute_object(|tc| {
            tc.mirna = vec![String::new(), "miRNA_stem".to_string()];
        });
        assert_eq!(past_the_end["mirna"], json!([null, "miRNA_stem"]));
        let none = attribute_object(|_| {});
        assert!(none.get("mirna").is_none());
    }

    /// The overlapping structural variants are one record-level list,
    /// `colocated_structural_variants` (JSON.pm:110 renames `sv`), in no
    /// consequence object (the factory builds those from an empty hash,
    /// JSON.pm:320), and absent when the record overlaps none; Ensembl VEP
    /// 116.3 writes `["nsv2770691","nsv3109934","nsv3163927","nsv3167562",
    /// "nsv533679"]` for `21 5028999 A G` of the `GRCh38-checksvs` corpus.
    #[test]
    fn colocated_structural_variants_is_a_record_level_list() {
        let mut v = make_annotated_variant();
        v.overlapping_svs = ["nsv2770691", "nsv3109934", "nsv533679"]
            .map(String::from)
            .to_vec();
        let formatter = fmt(FieldOptions {
            check_svs: true,
            ..Default::default()
        });
        let parsed = formatter.format_variant_value(&v);
        assert_eq!(
            parsed["colocated_structural_variants"],
            json!(["nsv2770691", "nsv3109934", "nsv533679"])
        );
        let consequence = &parsed["transcript_consequences"][0];
        for key in ["sv", "colocated_structural_variants"] {
            assert!(consequence.get(key).is_none(), "{key} in {consequence}");
        }
        v.overlapping_svs.clear();
        let parsed = formatter.format_variant_value(&v);
        assert!(parsed.get("colocated_structural_variants").is_none());
    }

    /// Under `--cell_type`, `cell_type` is an array of `<name>:<activity>`
    /// strings on a `regulatory_feature_consequences` entry, one element per
    /// requested epigenome the feature has an activity for whatever the name
    /// holds (a comma included), and the empty array on a
    /// `motif_feature_consequences` entry, which `JSON.pm:327` keeps, an
    /// arrayref being never `-`; a transcript entry has no such key, and
    /// without the flag no entry does.
    #[test]
    fn cell_type_is_an_array_and_empty_on_motif_rows() {
        use vep_core::consequence::{Activity, CellTypeActivity};
        let entry = |name: &str, activity| CellTypeActivity {
            name: name.into(),
            activity,
        };
        let mut v = make_annotated_variant();
        v.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.feature_id = Some("ENSR00000001".into());
            consequence.feature_biotype = Some("promoter".into());
            consequence.consequences = smallvec![Consequence::RegulatoryRegionVariant];
            consequence.impact = Impact::MODIFIER;
            consequence.feature_type = FeatureType::RegulatoryFeature;
            consequence.cell_type = Some(Box::new([
                entry("HUVEC", Activity::Repressed),
                entry("A549_(m,_58_y)", Activity::Active),
            ]));
            consequence
        });
        v.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.feature_id = Some("ENSM00000001".into());
            consequence.consequences = smallvec![Consequence::TfBindingSiteVariant];
            consequence.impact = Impact::MODIFIER;
            consequence.feature_type = FeatureType::MotifFeature;
            consequence.strand = 1;
            consequence.cell_type = Some(Box::new([]));
            consequence
        });
        let parsed = fmt(FieldOptions {
            regulatory: true,
            cell_type: true,
            ..Default::default()
        })
        .format_variant_value(&v);
        assert_eq!(
            parsed["regulatory_feature_consequences"][0]["cell_type"],
            json!(["HUVEC:REPRESSED", "A549_(m,_58_y):ACTIVE"])
        );
        assert_eq!(
            parsed["motif_feature_consequences"][0]["cell_type"],
            json!([])
        );
        assert!(parsed["transcript_consequences"][0]
            .get("cell_type")
            .is_none());
        for row in &mut v.transcript_consequences {
            row.cell_type = None;
        }
        let parsed = fmt(FieldOptions {
            regulatory: true,
            ..Default::default()
        })
        .format_variant_value(&v);
        for kind in [
            "transcript_consequences",
            "regulatory_feature_consequences",
            "motif_feature_consequences",
        ] {
            assert!(parsed[kind][0].get("cell_type").is_none(), "{kind}");
        }
    }

    /// A feature identifier and a transcription factor name are strings whatever
    /// their text, `007` and `1e3` included, while `motif_name` and the other
    /// members take the number conversion. Ensembl VEP exempts `transcript_id`
    /// alone of the feature keys (`JSON.pm:114-121`), so `numberify` (`:288`,
    /// `Utils.pm:231`, `:241`) writes an all-digit regulatory or motif
    /// identifier, and a factor named by digits, as a number (`007` as 7, `1e3`
    /// as 1000), the documented class `regulatory_identifier_numberified`.
    #[test]
    fn feature_identifiers_and_factor_names_stay_strings() {
        let regulatory = |id: &str, biotype: &str| {
            let mut consequence = TranscriptConsequence::default();
            consequence.feature_id = Some(id.into());
            consequence.feature_biotype = Some(biotype.into());
            consequence.consequences = smallvec![Consequence::RegulatoryRegionVariant];
            consequence.impact = Impact::MODIFIER;
            consequence.feature_type = FeatureType::RegulatoryFeature;
            consequence
        };
        let mut v = InputVariant::new(
            "21".into(),
            5_014_322,
            5_014_344,
            b"TGGGCCTGTGGGGGGCGCTTACC".to_vec(),
            b"-".to_vec(),
        );
        v.transcript_consequences = vec![
            regulatory("12345", "enhancer"),
            regulatory("007", "promoter_flanking_region"),
            regulatory("1e3", "open_chromatin_region"),
            regulatory("ENSR21_536SRH", "CTCF_binding_site"),
            {
                let mut row = TranscriptConsequence::default();
                row.feature_id = Some("424242".into());
                row.consequences = smallvec![Consequence::TfBindingSiteVariant];
                row.impact = Impact::MODIFIER;
                row.feature_type = FeatureType::MotifFeature;
                row.strand = 1;
                row.motif = Some(Box::new(vep_core::consequence::MotifFields::new(
                    Arc::from("12345"),
                    Some(7),
                    None,
                    None,
                    ["HNF4A", "123", "RXRA"]
                        .into_iter()
                        .map(Arc::from)
                        .collect(),
                )));
                row
            },
        ];
        v.most_severe_consequence = Some(Consequence::RegulatoryRegionVariant);
        let parsed = fmt(FieldOptions {
            regulatory: true,
            ..Default::default()
        })
        .format_variant_value(&v);
        let ids: Vec<&Value> = parsed["regulatory_feature_consequences"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| &row["regulatory_feature_id"])
            .collect();
        // The entries stand in identifier order, as the formats list a record's
        // regulatory rows.
        assert_eq!(
            ids,
            [
                &json!("007"),
                &json!("12345"),
                &json!("1e3"),
                &json!("ENSR21_536SRH")
            ]
        );
        let motif = &parsed["motif_feature_consequences"][0];
        assert_eq!(motif["motif_feature_id"], json!("424242"));
        assert_eq!(
            motif["transcription_factors"],
            json!(["HNF4A", "123", "RXRA"])
        );
        assert_eq!(motif["motif_name"], json!(12345));
        assert_eq!(motif["motif_pos"], json!(7));
    }

    #[test]
    fn no_header_for_json() {
        let mut f = fmt(FieldOptions::default());
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert!(buf.is_empty());
    }

    /// The consequence object of the annotated variant under `--sift`,
    /// `--polyphen` and `--domains`, with the given values on its consequence.
    fn prediction_object(edit: impl Fn(&mut TranscriptConsequence)) -> Value {
        let mut v = make_annotated_variant();
        edit(&mut v.transcript_consequences[0]);
        let options = FieldOptions {
            sift: true,
            polyphen: true,
            domains: true,
            ..Default::default()
        };
        fmt(options).format_variant_value(&v)["transcript_consequences"][0].clone()
    }

    /// `JSON.pm:353-360` splits each tool's value with
    /// `([a-z_]+)?\(?([\d\.]+)?\)?` and sets `<tool>_prediction` only for a
    /// prediction and `<tool>_score` only for a score: `--sift s` prints the bare
    /// score, so its object carries `sift_score` alone; `--sift p` carries
    /// `sift_prediction` alone; `b` carries both. The score is a number as
    /// `numberify` makes it, a whole score (`tolerated(1)`) an integer.
    #[test]
    fn score_only_mode_sets_score_key_only() {
        let score_only = prediction_object(|tc| {
            tc.sift = Some("0.01".to_string());
            tc.polyphen = Some("0.998".to_string());
        });
        assert_eq!(score_only["sift_score"], json!(0.01));
        assert_eq!(score_only["polyphen_score"], json!(0.998));
        assert!(score_only.get("sift_prediction").is_none(), "{score_only}");
        assert!(
            score_only.get("polyphen_prediction").is_none(),
            "{score_only}"
        );
        let prediction_only = prediction_object(|tc| {
            tc.sift = Some("deleterious_low_confidence".to_string());
            tc.polyphen = Some("probably_damaging".to_string());
        });
        assert_eq!(
            prediction_only["sift_prediction"],
            json!("deleterious_low_confidence")
        );
        assert_eq!(
            prediction_only["polyphen_prediction"],
            json!("probably_damaging")
        );
        assert!(
            prediction_only.get("sift_score").is_none(),
            "{prediction_only}"
        );
        assert!(
            prediction_only.get("polyphen_score").is_none(),
            "{prediction_only}"
        );
        let both = prediction_object(|tc| {
            tc.sift = Some("tolerated(1)".to_string());
            tc.polyphen = Some("benign(0.09)".to_string());
        });
        assert_eq!(both["sift_prediction"], json!("tolerated"));
        assert_eq!(both["sift_score"], json!(1));
        assert_eq!(both["polyphen_prediction"], json!("benign"));
        assert_eq!(both["polyphen_score"], json!(0.09));
        assert_eq!(
            parse_prediction_score("deleterious(0.01)"),
            (Some("deleterious"), Some("0.01"))
        );
        assert_eq!(parse_prediction_score("0"), (None, Some("0")));
    }

    /// `JSON.pm:363-371` re-splits each `DOMAINS` entry with `(\w+):(\w+)`, so
    /// the name stops at the first non-word character, and an entry whose `:`
    /// has no word character on both sides repeats the previous entry's split,
    /// `$1` and `$2` keeping the last successful match: the Ensembl VEP 116.3
    /// JSON for
    /// `DOMAINS=Gene3D:1.10.510.10,PANTHER:PTHR24416:SF8,AFDB-ENSP_mappings:AF-P60331-F1,Low_complexity_(Seg):seg,Pfam:PF00080`
    /// is `{Gene3D 1}` with the name a number, `numberify` (JSON.pm:288) having
    /// run over the record, `{PANTHER PTHR24416}`, `{ENSP_mappings AF}`,
    /// `{ENSP_mappings AF}` again and `{Pfam PF00080}`; a list opening with such
    /// an entry drops it (the GRCh38 chromosome 21 rows of ENST00000391624,
    /// `Low_complexity_(Seg):seg,AFDB-ENSP_mappings:AF-Q3LI67-F1`, carry one
    /// object). A capture of `0` fails `push ... if $1 && $2`, so
    /// `Gene3D:0.10.510.10,Low_complexity_(Seg):seg` is empty, and
    /// `Pfam:PF00080,Gene3D:0.10.510.10,Low_complexity_(Seg):seg` is
    /// `{Pfam PF00080}` alone, the `0` match having replaced the captures.
    #[test]
    fn domains_json_resplit_truncates_at_non_word() {
        let tc = prediction_object(|tc| {
            tc.domains = vec![
                ("Gene3D".to_string(), "1.10.510.10".to_string()),
                ("PANTHER".to_string(), "PTHR24416:SF8".to_string()),
                ("AFDB-ENSP_mappings".to_string(), "AF-P60331-F1".to_string()),
                ("Low_complexity_(Seg)".to_string(), "seg".to_string()),
                ("Pfam".to_string(), "PF00080".to_string()),
            ];
        });
        assert_eq!(
            tc["domains"],
            json!([
                {"db": "Gene3D", "name": 1},
                {"db": "PANTHER", "name": "PTHR24416"},
                {"db": "ENSP_mappings", "name": "AF"},
                {"db": "ENSP_mappings", "name": "AF"},
                {"db": "Pfam", "name": "PF00080"}
            ])
        );
        let opening = prediction_object(|tc| {
            tc.domains = vec![
                ("Low_complexity_(Seg)".to_string(), "seg".to_string()),
                ("AFDB-ENSP_mappings".to_string(), "AF-Q3LI67-F1".to_string()),
            ];
        });
        assert_eq!(
            opening["domains"],
            json!([{"db": "ENSP_mappings", "name": "AF"}])
        );
        // A part the text column drops (`grep {$_}`: empty or `0`) leaves an
        // entry with no `:`, which the pattern does not match.
        let dropped_part = prediction_object(|tc| {
            tc.domains = vec![
                ("Pfam".to_string(), "PF00080".to_string()),
                ("0".to_string(), "SSF46934".to_string()),
                (String::new(), "seg".to_string()),
            ];
        });
        assert_eq!(
            dropped_part["domains"],
            json!([
                {"db": "Pfam", "name": "PF00080"},
                {"db": "Pfam", "name": "PF00080"},
                {"db": "Pfam", "name": "PF00080"}
            ])
        );
        let zero_name = prediction_object(|tc| {
            tc.domains = vec![
                ("Gene3D".to_string(), "0.10.510.10".to_string()),
                ("Low_complexity_(Seg)".to_string(), "seg".to_string()),
            ];
        });
        assert_eq!(zero_name["domains"], json!([]));
        let zero_after_match = prediction_object(|tc| {
            tc.domains = vec![
                ("Pfam".to_string(), "PF00080".to_string()),
                ("Gene3D".to_string(), "0.10.510.10".to_string()),
                ("Low_complexity_(Seg)".to_string(), "seg".to_string()),
            ];
        });
        assert_eq!(
            zero_after_match["domains"],
            json!([{"db": "Pfam", "name": "PF00080"}])
        );
        assert_eq!(
            word_colon_word("a:b:c"),
            Some(("a".to_string(), "b".to_string()))
        );
        assert_eq!(word_colon_word(":x"), None);
        assert_eq!(word_colon_word("x:"), None);
    }

    /// A row carrying the AlphaMissense plugin's fields writes them under
    /// `alphamissense` with their names as keys (AlphaMissense.pm:281-282,
    /// JSON.pm:345), the score as a number and the class as text; no flat `am_*`
    /// key is written, and a row without the fields has no `alphamissense` key.
    #[test]
    fn alphamissense_json_nests_the_selected_columns() {
        let plugin_values = [
            ("am_pathogenicity", "0.9123"),
            ("am_class", "likely_pathogenic"),
            ("am_transcript_id", "ENST00000335137.4"),
        ];
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_values
                .iter()
                .map(|(key, _)| key.to_string())
                .collect(),
            Some("GRCh38".into()),
        );
        let mut v = make_annotated_variant();
        for (key, value) in plugin_values {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), value.to_string());
        }
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["alphamissense"],
            json!({
                "am_pathogenicity": 0.9123,
                "am_class": "likely_pathogenic",
                "am_transcript_id": "ENST00000335137.4",
            })
        );
        let flat: Vec<&String> = row
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| key.starts_with("am_"))
            .collect();
        assert!(flat.is_empty(), "{flat:?}");
        let bare = formatter.format_variant_value(&make_annotated_variant());
        assert!(bare["transcript_consequences"][0]
            .get("alphamissense")
            .is_none());
    }

    /// A row carrying the SpliceVault plugin's fields writes them once under
    /// `splicevault` with their names as keys (SpliceVault.pm:191, JSON.pm:345),
    /// the counts and the delta as numbers, and the events as an object keyed
    /// by rank whose impact is a number when it reads as one (`+48` is 48 to
    /// `numberify`), whose percent stays text, and whose missing parts are
    /// `null` (SpliceVault.pm:180-191); no flat `splicevault_*` key is written,
    /// and a row without the fields has no `splicevault` key. An empty event of
    /// the plugin's list (`a||b`, which the text formats keep as `a,,b`) is keyed
    /// `""` with four `null` parts, and the empty list is `{}`.
    /// The row's cells are made up.
    #[test]
    fn splicevault_json_nests_rank_keyed_events() {
        let plugin_values = [
            ("SpliceVault_site_sample_count", "16054"),
            ("SpliceVault_site_max_depth", "17"),
            ("SpliceVault_SpliceAI_delta", "0.36"),
            ("SpliceVault_out_of_frame_events", "Frameshift:2/2"),
            (
                "SpliceVault_top_events",
                "Top1:CD:-35:31%:Frameshift|Top2:CD:+48:0.006%:Frameshift|\
                 Top3:ES:5-6:15_840%:inFrame|Top4:CA:+7",
            ),
            ("SpliceVault_site_pos", "chr1:925801"),
            ("SpliceVault_site_type", "Donor_loss"),
        ];
        let plugin_fields: Vec<String> = plugin_values
            .iter()
            .map(|(key, _)| key.to_string())
            .collect();
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_fields,
            Some("GRCh38".into()),
        );
        let mut v = make_annotated_variant();
        for (key, value) in plugin_values {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), value.to_string());
        }
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["splicevault"],
            json!({
                "SpliceVault_site_sample_count": 16054,
                "SpliceVault_site_max_depth": 17,
                "SpliceVault_SpliceAI_delta": 0.36,
                "SpliceVault_out_of_frame_events": "Frameshift:2/2",
                "SpliceVault_top_events": {
                    "Top1": {"site_type": "CD", "transcript_impact": -35,
                             "percent_of_supporting_samples": "31%", "frame": "Frameshift"},
                    "Top2": {"site_type": "CD", "transcript_impact": 48,
                             "percent_of_supporting_samples": "0.006%", "frame": "Frameshift"},
                    "Top3": {"site_type": "ES", "transcript_impact": "5-6",
                             "percent_of_supporting_samples": "15_840%", "frame": "inFrame"},
                    "Top4": {"site_type": "CA", "transcript_impact": 7,
                             "percent_of_supporting_samples": null, "frame": null},
                },
                "SpliceVault_site_pos": "chr1:925801",
                "SpliceVault_site_type": "Donor_loss",
            })
        );
        let flat: Vec<&String> = row
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| key.starts_with("splicevault_"))
            .collect();
        assert!(flat.is_empty(), "{flat:?}");
        let bare = formatter.format_variant_value(&make_annotated_variant());
        assert!(bare["transcript_consequences"][0]
            .get("splicevault")
            .is_none());
        assert_eq!(
            splicevault_events_object(""),
            json!({}),
            "no event is an empty object, as `split` of an empty list leaves it"
        );
        let unranked = json!({"site_type": null, "transcript_impact": null,
                              "percent_of_supporting_samples": null, "frame": null});
        assert_eq!(
            splicevault_events_object("Top1:CD:-35:31%:Frameshift||Top2:ES:5-6:31%:inFrame"),
            json!({
                "Top1": {"site_type": "CD", "transcript_impact": -35,
                         "percent_of_supporting_samples": "31%", "frame": "Frameshift"},
                "": unranked,
                "Top2": {"site_type": "ES", "transcript_impact": "5-6",
                         "percent_of_supporting_samples": "31%", "frame": "inFrame"},
            }),
            "an interior empty event is keyed by its undefined rank"
        );
        assert_eq!(
            splicevault_events_object("|Top1:CD:-35:31%:Frameshift"),
            json!({
                "": unranked,
                "Top1": {"site_type": "CD", "transcript_impact": -35,
                         "percent_of_supporting_samples": "31%", "frame": "Frameshift"},
            }),
            "a leading empty event is keyed the same way"
        );
    }

    /// `split(/:/, $event)` (SpliceVault.pm:184) splits an event on `:` alone,
    /// so a `,` anywhere in an event is text: inside the percent (`1,5%`), in
    /// the frame, or in the rank, where it is part of the key. Each such event
    /// is one event, and no event is made of the text after the `,`.
    /// The events are made up.
    #[test]
    fn splicevault_event_holding_a_comma_is_one_event() {
        assert_eq!(
            splicevault_events_object("Top1:CD:-330:73_314%:inFrame|Top2:CA:+112:1,5%:inFrame"),
            json!({
                "Top1": {"site_type": "CD", "transcript_impact": -330,
                         "percent_of_supporting_samples": "73_314%", "frame": "inFrame"},
                "Top2": {"site_type": "CA", "transcript_impact": 112,
                         "percent_of_supporting_samples": "1,5%", "frame": "inFrame"},
            })
        );
        assert_eq!(
            splicevault_events_object("Top1:CD:-330:73%:in,Frame|To,p2:CA:1,5:7%:inFrame"),
            json!({
                "Top1": {"site_type": "CD", "transcript_impact": -330,
                         "percent_of_supporting_samples": "73%", "frame": "in,Frame"},
                "To,p2": {"site_type": "CA", "transcript_impact": "1,5",
                          "percent_of_supporting_samples": "7%", "frame": "inFrame"},
            }),
            "a `,` in the frame or the rank is text too, and `1,5` is no number"
        );
    }

    /// Inside the `splicevault` object a cell the table row lacks is `null`
    /// (the twelve-variable split of SpliceVault.pm:201 leaves it undefined and
    /// the encoder writes the nested hash as it stands; JSON.pm:327 reads the
    /// top-level keys alone) and a cell the row has and leaves empty, which the
    /// plugin names in `plugin_empty_strings`, is `""`; neither writes a flat
    /// `splicevault_*` key, and an empty events cell is `{}` either way.
    /// The row's cells are made up.
    #[test]
    fn splicevault_json_writes_null_for_a_missing_cell_and_the_empty_string_for_an_empty_one() {
        let plugin_values = [
            ("SpliceVault_site_sample_count", "16054"),
            ("SpliceVault_site_max_depth", ""),
            ("SpliceVault_SpliceAI_delta", ""),
            ("SpliceVault_out_of_frame_events", "Frameshift:2/2"),
            ("SpliceVault_top_events", ""),
            ("SpliceVault_site_pos", "chr1:925801"),
            ("SpliceVault_site_type", "Donor_loss"),
        ];
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_values
                .iter()
                .map(|(key, _)| key.to_string())
                .collect(),
            Some("GRCh38".into()),
        );
        let mut v = make_annotated_variant();
        for (key, value) in plugin_values {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), value.to_string());
        }
        v.transcript_consequences[0]
            .plugin_empty_strings
            .push("SpliceVault_SpliceAI_delta".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["splicevault"],
            json!({
                "SpliceVault_site_sample_count": 16054,
                "SpliceVault_site_max_depth": null,
                "SpliceVault_SpliceAI_delta": "",
                "SpliceVault_out_of_frame_events": "Frameshift:2/2",
                "SpliceVault_top_events": {},
                "SpliceVault_site_pos": "chr1:925801",
                "SpliceVault_site_type": "Donor_loss",
            })
        );
        let flat: Vec<&String> = row
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| key.starts_with("splicevault_"))
            .collect();
        assert!(flat.is_empty(), "{flat:?}");
    }

    /// `split(/:/, $event)` into five variables (SpliceVault.pm:184) runs with
    /// a limit of six, which keeps a trailing empty part: an event ending in
    /// `:` has the empty string for its `frame`, one ending in `::` for its
    /// percent as well, while a part the event lacks stays `null` and a sixth
    /// part and anything after it are dropped.
    /// The events are made up.
    #[test]
    fn splicevault_event_trailing_empty_parts_are_empty_strings() {
        assert_eq!(
            splicevault_events_object("Top1:CA:-480:7%:|Top2:CA:-76:39%:Frameshift"),
            json!({
                "Top1": {"site_type": "CA", "transcript_impact": -480,
                         "percent_of_supporting_samples": "7%", "frame": ""},
                "Top2": {"site_type": "CA", "transcript_impact": -76,
                         "percent_of_supporting_samples": "39%", "frame": "Frameshift"},
            })
        );
        assert_eq!(
            splicevault_events_object("Top1:CD:-330::"),
            json!({
                "Top1": {"site_type": "CD", "transcript_impact": -330,
                         "percent_of_supporting_samples": "", "frame": ""},
            })
        );
        assert_eq!(
            splicevault_events_object("Top1::::"),
            json!({
                "Top1": {"site_type": "", "transcript_impact": "",
                         "percent_of_supporting_samples": "", "frame": ""},
            })
        );
        let five_parts = json!({
            "Top1": {"site_type": "CD", "transcript_impact": -330,
                     "percent_of_supporting_samples": "73%", "frame": "inFrame"},
        });
        for event in [
            "Top1:CD:-330:73%:inFrame:",
            "Top1:CD:-330:73%:inFrame::",
            "Top1:CD:-330:73%:inFrame:x:y",
        ] {
            assert_eq!(splicevault_events_object(event), five_parts, "{event}");
        }
        assert_eq!(
            splicevault_events_object("Top3:ES:7"),
            json!({
                "Top3": {"site_type": "ES", "transcript_impact": 7,
                         "percent_of_supporting_samples": null, "frame": null},
            })
        );
    }

    /// A row carrying the SpliceAI plugin's nine fields and its flag writes them
    /// once under `spliceai` with the `SpliceAI_pred_` prefix dropped
    /// (SpliceAI.pm:275-288, 316, JSON.pm:345): the scores and positions as
    /// numbers (`0.00` is 0, `-22` is -22), `SYMBOL` a string whatever it looks
    /// like (JSON.pm:120) and the flag under its own name; no flat
    /// `spliceai_pred_*` key is written, a row without the fields has no
    /// `spliceai` key, and a row without the flag has no `SpliceAI_cutoff` entry.
    /// The row's cells are made up.
    #[test]
    fn spliceai_json_nests_the_parts_under_the_plugin_name() {
        let plugin_values = [
            ("SpliceAI_pred_SYMBOL", "1060P11.3"),
            ("SpliceAI_pred_DS_AG", "0.00"),
            ("SpliceAI_pred_DS_AL", "0.01"),
            ("SpliceAI_pred_DS_DG", "0.00"),
            ("SpliceAI_pred_DS_DL", "0.42"),
            ("SpliceAI_pred_DP_AG", "2"),
            ("SpliceAI_pred_DP_AL", "37"),
            ("SpliceAI_pred_DP_DG", "-22"),
            ("SpliceAI_pred_DP_DL", "37"),
            ("SpliceAI_cutoff", "FAIL"),
        ];
        let plugin_fields: Vec<String> = plugin_values
            .iter()
            .map(|(key, _)| key.to_string())
            .collect();
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_fields,
            Some("GRCh38".into()),
        );
        let row_with = |values: &[(&str, &str)]| {
            let mut v = make_annotated_variant();
            for (key, value) in values {
                v.transcript_consequences[0]
                    .plugin_data
                    .insert(key.to_string(), value.to_string());
            }
            formatter.format_variant_value(&v)["transcript_consequences"][0].clone()
        };
        let row = row_with(&plugin_values);
        assert_eq!(
            row["spliceai"],
            json!({
                "SYMBOL": "1060P11.3",
                "DS_AG": 0, "DS_AL": 0.01, "DS_DG": 0, "DS_DL": 0.42,
                "DP_AG": 2, "DP_AL": 37, "DP_DG": -22, "DP_DL": 37,
                "SpliceAI_cutoff": "FAIL",
            })
        );
        let flat: Vec<&String> = row
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| key.starts_with("spliceai_"))
            .collect();
        assert!(flat.is_empty(), "{flat:?}");
        let bare = row_with(&[]);
        assert!(bare.get("spliceai").is_none());
        let without_flag = row_with(&plugin_values[..9]);
        assert_eq!(without_flag["spliceai"].as_object().unwrap().len(), 9);
        assert!(without_flag["spliceai"].get("SpliceAI_cutoff").is_none());
    }

    /// A row whose gene several entries of one SpliceAI score record name
    /// carries each member's values joined by commas: the member is that text,
    /// a string, where a member of one value is a number, and the flag stays
    /// one verdict. The row's cells are made up.
    #[test]
    fn spliceai_json_writes_a_member_of_several_values_as_its_text() {
        let plugin_values = [
            ("SpliceAI_pred_SYMBOL", "GENE1,GENE1"),
            ("SpliceAI_pred_DS_AG", "0.10,0.11"),
            ("SpliceAI_pred_DS_AL", "0.20,0.21"),
            ("SpliceAI_pred_DS_DG", "0.30,0.31"),
            ("SpliceAI_pred_DS_DL", "0.40,0.41"),
            ("SpliceAI_pred_DP_AG", "1,-1"),
            ("SpliceAI_pred_DP_AL", "2,-2"),
            ("SpliceAI_pred_DP_DG", "3,-3"),
            ("SpliceAI_pred_DP_DL", "4,-4"),
            ("SpliceAI_cutoff", "PASS"),
        ];
        let plugin_fields: Vec<String> = plugin_values
            .iter()
            .map(|(key, _)| key.to_string())
            .collect();
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_fields,
            Some("GRCh38".into()),
        );
        let mut v = make_annotated_variant();
        for (key, value) in plugin_values {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), value.to_string());
        }
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["spliceai"],
            json!({
                "SYMBOL": "GENE1,GENE1",
                "DS_AG": "0.10,0.11", "DS_AL": "0.20,0.21",
                "DS_DG": "0.30,0.31", "DS_DL": "0.40,0.41",
                "DP_AG": "1,-1", "DP_AL": "2,-2", "DP_DG": "3,-3", "DP_DL": "4,-4",
                "SpliceAI_cutoff": "PASS",
            })
        );
    }

    /// Inside the `spliceai` and `alphamissense` objects an undefined value, the
    /// part that a truncated SpliceAI entry lacks or the cell that a short
    /// AlphaMissense row lacks, is `null`, the empty string a row's empty cell
    /// gives is `""` and a `-` cell is `"-"`, as the encoder writes the plugins'
    /// nested hashes as they stand (SpliceAI.pm:277-287 reads nine parts into
    /// nine keys, AlphaMissense.pm:279 every selected column);
    /// `OutputFactory/JSON.pm:327` drops an undefined or `-` top-level key
    /// alone. No flat key is written for either.
    /// The row's cells are made up.
    #[test]
    fn nested_plugin_objects_write_an_undefined_value_as_null_and_an_empty_cell_as_a_string() {
        let plugin_fields: Vec<String> = [
            "SpliceAI_pred_SYMBOL",
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
            "SpliceAI_pred_DP_AG",
            "SpliceAI_pred_DP_AL",
            "SpliceAI_pred_DP_DG",
            "SpliceAI_pred_DP_DL",
            "am_pathogenicity",
            "am_class",
        ]
        .iter()
        .map(|key| key.to_string())
        .collect();
        let formatter = JsonOutputFormatter::new(
            FieldOptions::default(),
            plugin_fields.clone(),
            Some("GRCh38".into()),
        );
        let mut v = make_annotated_variant();
        let tc = &mut v.transcript_consequences[0];
        for key in &plugin_fields {
            tc.plugin_data.insert(key.clone(), String::new());
        }
        tc.plugin_data
            .insert("SpliceAI_pred_SYMBOL".to_string(), "BTG3-AS1".to_string());
        tc.plugin_data
            .insert("SpliceAI_pred_DS_AG".to_string(), "0.50".to_string());
        tc.plugin_data
            .insert("SpliceAI_pred_DS_DG".to_string(), "0.49".to_string());
        tc.plugin_data
            .insert("am_class".to_string(), "likely_benign".to_string());
        tc.plugin_empty_strings = vec![
            "SpliceAI_pred_DS_AL".to_string(),
            "am_pathogenicity".to_string(),
        ];
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["spliceai"],
            json!({
                "SYMBOL": "BTG3-AS1",
                "DS_AG": 0.5, "DS_AL": "", "DS_DG": 0.49, "DS_DL": null,
                "DP_AG": null, "DP_AL": null, "DP_DG": null, "DP_DL": null,
            })
        );
        assert_eq!(
            row["alphamissense"],
            json!({"am_pathogenicity": "", "am_class": "likely_benign"})
        );
        v.transcript_consequences[0]
            .plugin_data
            .insert("am_class".to_string(), String::new());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["alphamissense"],
            json!({"am_pathogenicity": "", "am_class": null})
        );
        v.transcript_consequences[0]
            .plugin_data
            .insert("am_class".to_string(), "-".to_string());
        let row = formatter.format_variant_value(&v)["transcript_consequences"][0].clone();
        assert_eq!(
            row["alphamissense"],
            json!({"am_pathogenicity": "", "am_class": "-"})
        );
        let flat: Vec<&String> = row
            .as_object()
            .unwrap()
            .keys()
            .filter(|key| key.starts_with("spliceai_") || key.starts_with("am_"))
            .collect();
        assert!(flat.is_empty(), "{flat:?}");
    }

    /// A plugin row carrying the UTRAnnotator keys: `5UTR_annotation` nests one
    /// object per `,`-separated group under `1`, `2`, ... with its numbers as
    /// numbers and `NA`, `-` and the empty strength as strings; the flag and the
    /// counts take the generic lowercased keys; a `-` value is left out as every
    /// other `-` is.
    #[test]
    fn utr_annotation_nests_one_object_per_group() {
        let plugin_object = |annotation: &str| {
            let mut v = make_annotated_variant();
            let tc = &mut v.transcript_consequences[0];
            tc.plugin_data.insert(
                "5UTR_consequence".to_string(),
                "5_prime_UTR_uORF_frameshift_variant".to_string(),
            );
            tc.plugin_data
                .insert("5UTR_annotation".to_string(), annotation.to_string());
            tc.plugin_data
                .insert("Existing_uORFs".to_string(), "2".to_string());
            let formatter = JsonOutputFormatter::new(
                FieldOptions::default(),
                vec![
                    "5UTR_annotation".to_string(),
                    "5UTR_consequence".to_string(),
                    "Existing_uORFs".to_string(),
                ],
                Some("GRCh38".into()),
            );
            formatter.format_variant_value(&v)["transcript_consequences"][0].clone()
        };
        let two = plugin_object(
            "Evidence=NA:KozakContext=CCCATGC:KozakStrength=Weak:alt_type=InFrame_oORF:alt_type_length=111:ref_StartDistanceToCDS=20:ref_type=OutOfFrame_oORF:ref_type_length=NA,\
             Evidence=False:KozakContext=-:KozakStrength=:alt_type=uORF:alt_type_length=18:ref_StartDistanceToCDS=12:ref_type=uORF:ref_type_length=0",
        );
        assert_eq!(
            two["5utr_annotation"],
            json!({
                "1": {"Evidence": "NA", "KozakContext": "CCCATGC", "KozakStrength": "Weak",
                      "alt_type": "InFrame_oORF", "alt_type_length": 111,
                      "ref_StartDistanceToCDS": 20, "ref_type": "OutOfFrame_oORF",
                      "ref_type_length": "NA"},
                "2": {"Evidence": "False", "KozakContext": "-", "KozakStrength": "",
                      "alt_type": "uORF", "alt_type_length": 18,
                      "ref_StartDistanceToCDS": 12, "ref_type": "uORF",
                      "ref_type_length": 0}
            })
        );
        assert_eq!(
            two["5utr_consequence"],
            "5_prime_UTR_uORF_frameshift_variant"
        );
        assert_eq!(two["existing_uorfs"], 2);
        let none = plugin_object("-");
        assert!(none.get("5utr_annotation").is_none(), "{none}");
    }

    /// The JSON factory moves the custom annotations to the record before any
    /// row hash is built (OutputFactory/JSON.pm:249-258), so a `--custom`
    /// source named `SOURCE` lists its names under `custom_annotations` and no
    /// entry carries a `source` from it; an entry's `source` is its
    /// transcript's source cache alone (OutputFactory.pm:1592).
    #[test]
    fn a_custom_source_named_source_writes_no_source_on_an_entry() {
        use crate::output::fields::CustomSourceInfo;
        use vep_core::consequence::SourceCache;
        use vep_core::custom_annotation::{CustomAnnotations, CustomRecord};
        let annotations = vec![CustomAnnotations::new(
            Arc::from("SOURCE"),
            vec![CustomRecord::new(
                "exact_1".to_string(),
                None,
                indexmap::IndexMap::new(),
                None,
            )],
            true,
            None,
        )];
        let options = FieldOptions {
            custom_sources: vec![CustomSourceInfo {
                file: "/data/regions.bed.gz".to_string(),
                short_name: "SOURCE".to_string(),
                fields: Vec::new(),
                filter_field: false,
                summary_stats: Vec::new(),
                reciprocal: false,
                field_descriptions: indexmap::IndexMap::new(),
            }],
            ..Default::default()
        };
        let mut plain = make_annotated_variant();
        plain.custom_annotations = annotations.clone();
        let value = fmt(options.clone()).format_variant_value(&plain);
        assert_eq!(
            value["custom_annotations"],
            json!({"SOURCE": [{"name": "exact_1"}]})
        );
        assert!(
            value["transcript_consequences"][0].get("source").is_none(),
            "{value}"
        );

        let mut merged = make_annotated_variant();
        merged.custom_annotations = annotations.clone();
        merged.transcript_consequences[0].transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.source_cache = Some(SourceCache::RefSeq);
            transcript
        }));
        let value = fmt(options.clone()).format_variant_value(&merged);
        assert_eq!(value["transcript_consequences"][0]["source"], "RefSeq");

        let mut intergenic = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        intergenic.custom_annotations = annotations;
        let value = fmt(options).format_variant_value(&intergenic);
        assert!(
            value["intergenic_consequences"][0].get("source").is_none(),
            "{value}"
        );
    }

    /// The same holds for a `--custom` source whose key spells any other flag
    /// field: the text formats print the source's values under the key where
    /// the field sets none, and here they stand under `custom_annotations`
    /// alone (OutputFactory/JSON.pm:249-258). A source named `SYMBOL` leaves
    /// an entry's `gene_symbol` its transcript's symbol or absent under
    /// `--symbol`, and one named `MINIMISED` writes no `minimised` on a record
    /// `--minimal` leaves as it is.
    #[test]
    fn a_custom_key_named_as_a_flag_field_writes_no_value_on_an_entry() {
        use crate::output::fields::{self, CustomSourceInfo};
        use vep_core::custom_annotation::{CustomAnnotations, CustomRecord};
        let annotations = |short_name: &str| {
            vec![CustomAnnotations::new(
                Arc::from(short_name),
                vec![CustomRecord::new(
                    "exact_1".to_string(),
                    None,
                    indexmap::IndexMap::new(),
                    None,
                )],
                true,
                None,
            )]
        };
        let source = |short_name: &str| CustomSourceInfo {
            file: "/data/regions.bed.gz".to_string(),
            short_name: short_name.to_string(),
            fields: Vec::new(),
            filter_field: false,
            summary_stats: Vec::new(),
            reciprocal: false,
            field_descriptions: indexmap::IndexMap::new(),
        };

        let options = FieldOptions {
            symbol: true,
            custom_sources: vec![source("SYMBOL")],
            ..Default::default()
        };
        let mut named = make_annotated_variant();
        named.custom_annotations = annotations("SYMBOL");
        let value = fmt(options.clone()).format_variant_value(&named);
        assert_eq!(
            value["custom_annotations"],
            json!({"SYMBOL": [{"name": "exact_1"}]})
        );
        assert_eq!(value["transcript_consequences"][0]["gene_symbol"], "MRPL39");

        let mut intergenic = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        intergenic.custom_annotations = annotations("SYMBOL");
        assert_eq!(
            fields::field_value("SYMBOL", &intergenic, None, &options, ","),
            "exact_1",
            "the text formats print the source's names under the key"
        );
        let value = fmt(options).format_variant_value(&intergenic);
        assert!(
            value["intergenic_consequences"][0]
                .get("gene_symbol")
                .is_none(),
            "{value}"
        );
        assert_eq!(
            value["custom_annotations"],
            json!({"SYMBOL": [{"name": "exact_1"}]})
        );

        let options = FieldOptions {
            minimal: true,
            custom_sources: vec![source("MINIMISED")],
            ..Default::default()
        };
        intergenic.custom_annotations = annotations("MINIMISED");
        assert_eq!(
            fields::field_value("MINIMISED", &intergenic, None, &options, ","),
            "exact_1"
        );
        let value = fmt(options).format_variant_value(&intergenic);
        assert!(value.get("minimised").is_none(), "{value}");

        // The empty list a regulatory row is assigned under `--cell_type`
        // (OutputFactory.pm:1888) is the field's own value and stays on the
        // entry beside a source named `CELL_TYPE`.
        let options = FieldOptions {
            regulatory: true,
            cell_type: true,
            custom_sources: vec![source("CELL_TYPE")],
            ..Default::default()
        };
        let mut regulatory = intergenic.clone();
        regulatory.custom_annotations = annotations("CELL_TYPE");
        regulatory.transcript_consequences = vec![{
            let mut consequence = TranscriptConsequence::default();
            consequence.feature_id = Some("ENSR00000001".into());
            consequence.consequences = smallvec![Consequence::RegulatoryRegionVariant];
            consequence.impact = Impact::MODIFIER;
            consequence.feature_type = FeatureType::RegulatoryFeature;
            consequence.cell_type = Some(Box::new([]));
            consequence
        }];
        let value = fmt(options).format_variant_value(&regulatory);
        assert_eq!(
            value["regulatory_feature_consequences"][0]["cell_type"],
            json!([]),
            "{value}"
        );

        // A flag field a plugin declares and leaves unwritten on the row is
        // read the same way: beside a source named `TSL`, the entry of a
        // transcript without a support level carries no `tsl`, with `--tsl`
        // as without it.
        named.custom_annotations = annotations("TSL");
        for tsl in [true, false] {
            let options = FieldOptions {
                tsl,
                custom_sources: vec![source("TSL")],
                ..Default::default()
            };
            let declaring =
                JsonOutputFormatter::new(options, vec!["TSL".to_string()], Some("GRCh37".into()));
            let value = declaring.format_variant_value(&named);
            assert!(
                value["transcript_consequences"][0].get("tsl").is_none(),
                "--tsl {tsl}: {value}"
            );
            assert_eq!(
                value["custom_annotations"],
                json!({"TSL": [{"name": "exact_1"}]})
            );
        }

        // A frequency column a plugin and a source both declare: with the
        // flag off the text formats print the source's names under it on a
        // row the plugin left unwritten, and the entry carries no such key.
        let mut options = FieldOptions {
            custom_sources: vec![source("gnomADg_AF")],
            ..Default::default()
        };
        options.set_plugin_fields(["gnomADg_AF".to_string()]);
        intergenic.custom_annotations = annotations("gnomADg_AF");
        assert_eq!(
            fields::field_value("gnomADg_AF", &intergenic, None, &options, ","),
            "exact_1"
        );
        let declaring = JsonOutputFormatter::new(
            options,
            vec!["gnomADg_AF".to_string()],
            Some("GRCh37".into()),
        );
        let value = declaring.format_variant_value(&intergenic);
        let entry = value["intergenic_consequences"][0]
            .as_object()
            .expect("the entry");
        assert!(
            !entry
                .keys()
                .any(|key| key.eq_ignore_ascii_case("gnomADg_AF")),
            "{value}"
        );
    }

    /// A flag field Ensembl VEP sets as a plain value before a `--custom`
    /// source's values are pushed under the same key holds, in the text
    /// formats, its value and the pushed values after it. The JSON factory
    /// pushes none (OutputFactory/JSON.pm:249-258), so a record's
    /// `variant_class` and an entry's `impact` and `allele_num` are the
    /// field's own, and the source's record stands under `custom_annotations`.
    #[test]
    fn a_custom_key_named_as_a_plain_value_field_writes_the_fields_value_alone_on_an_entry() {
        use crate::output::fields::{self, CustomSourceInfo};
        use vep_core::custom_annotation::{CustomAnnotations, CustomRecord};
        let annotations = |short_name: &str, allele: Option<&[u8]>| {
            vec![CustomAnnotations::new(
                Arc::from(short_name),
                vec![CustomRecord::new(
                    "exact_1".to_string(),
                    allele.map(<[u8]>::to_vec),
                    indexmap::IndexMap::new(),
                    None,
                )],
                true,
                None,
            )]
        };
        let source = |short_name: &str| CustomSourceInfo {
            file: "/data/sites.vcf.gz".to_string(),
            short_name: short_name.to_string(),
            fields: Vec::new(),
            filter_field: false,
            summary_stats: Vec::new(),
            reciprocal: false,
            field_descriptions: indexmap::IndexMap::new(),
        };
        let mut intergenic = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);

        // A value of the record's hash, beside a record of the whole input
        // record.
        let options = FieldOptions {
            variant_class: true,
            custom_sources: vec![source("VARIANT_CLASS")],
            ..Default::default()
        };
        intergenic.custom_annotations = annotations("VARIANT_CLASS", None);
        assert_eq!(
            fields::field_value("VARIANT_CLASS", &intergenic, None, &options, ","),
            "SNV,exact_1",
            "the text formats list the source's name after the class"
        );
        assert_eq!(
            fields::field_value_without_custom_pushes(
                "VARIANT_CLASS",
                &intergenic,
                None,
                &options,
                ","
            ),
            "SNV"
        );
        let value = fmt(options).format_variant_value(&intergenic);
        assert_eq!(value["variant_class"], "SNV", "{value}");
        assert_eq!(
            value["custom_annotations"],
            json!({"VARIANT_CLASS": [{"name": "exact_1"}]})
        );

        // Two values a row assigns, beside a record of the row's allele.
        for (key, json_key, text, own) in [
            ("IMPACT", "impact", "MODIFIER,exact_1", json!("MODIFIER")),
            ("ALLELE_NUM", "allele_num", "1,exact_1", json!(1)),
        ] {
            let options = FieldOptions {
                allele_number: true,
                custom_sources: vec![source(key)],
                ..Default::default()
            };
            intergenic.custom_annotations = annotations(key, Some(b"G"));
            assert_eq!(
                fields::field_value(key, &intergenic, None, &options, ","),
                text,
                "{key}"
            );
            let value = fmt(options).format_variant_value(&intergenic);
            assert_eq!(
                value["intergenic_consequences"][0][json_key], own,
                "{key}: {value}"
            );
            assert_eq!(
                value["custom_annotations"][key],
                json!([{"allele": "G", "name": "exact_1"}]),
                "{key}"
            );
        }
    }

    /// A field a plugin declares under the name of a co-located variant field
    /// (a frequency column, `CLIN_SIG`) or of a field of the record's own
    /// hash (`VARIANT_CLASS`) is on an entry only as the plugin's own. The
    /// JSON factory builds a row's hash from an empty one, the record's keys
    /// staying on the record (OutputFactory/JSON.pm:261-262, 320), sets no
    /// co-located key on it (:493-515), and a plugin's keys are written over
    /// the hash of a row the plugin ran on (OutputFactory.pm:2228-2241). So
    /// under the field's flag the entry of a row the plugin left unwritten
    /// carries no such key, where the text formats print the record's or the
    /// co-located value in the plugin's column; the record's own key and
    /// `colocated_variants` hold the data either way.
    #[test]
    fn a_plugin_field_no_json_row_holds_is_on_an_entry_only_where_the_plugin_wrote_it() {
        use crate::output::fields;
        let declared = ["gnomADg_AF", "ASN_AF", "CLIN_SIG", "VARIANT_CLASS"];
        let mut v = make_annotated_variant();
        let mut unwritten = v.transcript_consequences[0].clone();
        unwritten.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = "ENST00000307301".into();
            transcript
        }));
        for (key, value) in [
            ("gnomADg_AF", "0.61"),
            ("ASN_AF", "0.5"),
            ("CLIN_SIG", "from_the_table"),
            ("VARIANT_CLASS", "from_the_table"),
        ] {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), value.to_string());
        }
        v.transcript_consequences.push(unwritten);
        v.colocated_variants = vec![{
            let mut colocated = vep_core::variant::ColocatedVariant::default();
            colocated.chr = Some("21".into());
            colocated.id = "rs123".into();
            colocated.start = 25000100;
            colocated.end = 25000100;
            colocated.allele_string = Some("A/G".into());
            colocated.strand = 1;
            colocated.clin_sig = Some("uncertain_significance".into());
            colocated.matched_alleles = vec![vep_core::variant::MatchedAllele {
                input_allele: "G".to_string(),
                input_index: 0,
                record_allele: "G".to_string(),
                record_index: 0,
            }];
            colocated.frequency_strings = [("gnomADg", "G:0.0004729"), ("ASN", "G:0.25")]
                .into_iter()
                .map(|(column, value)| (column.to_string(), value.to_string()))
                .collect();
            colocated
        }];
        v.existing_variation = vec!["rs123".into()];

        let flags_on = FieldOptions {
            check_existing: true,
            af_1kg: true,
            af_gnomadg: true,
            variant_class: true,
            ..Default::default()
        };
        let flags_off = FieldOptions::default();
        for (options, flagged) in [(flags_on, true), (flags_off, false)] {
            let formatter = JsonOutputFormatter::new(
                options,
                declared.map(str::to_string).to_vec(),
                Some("GRCh38".into()),
            );
            let value = formatter.format_variant_value(&v);
            let entry_of = |transcript: &str| {
                value["transcript_consequences"]
                    .as_array()
                    .expect("the transcript entries")
                    .iter()
                    .find(|entry| entry["transcript_id"] == transcript)
                    .unwrap_or_else(|| panic!("no entry of {transcript}: {value}"))
            };
            let written = entry_of("ENST00000352957");
            assert_eq!(written["gnomadg_af"], json!(0.61), "{value}");
            assert_eq!(written["asn_af"], json!(0.5), "{value}");
            assert_eq!(written["clin_sig"], "from_the_table", "{value}");
            assert_eq!(written["variant_class"], "from_the_table", "{value}");
            let entry = entry_of("ENST00000307301")
                .as_object()
                .expect("the entry the plugin left unwritten");
            for key in ["gnomadg_af", "asn_af", "clin_sig", "variant_class"] {
                assert!(
                    !entry.contains_key(key),
                    "{key} on the entry the plugin left unwritten, flags on: {flagged}: {value}"
                );
            }
            assert_eq!(
                value["colocated_variants"][0]["frequencies"],
                json!({"G": {"asn": 0.25, "gnomadg": 0.0004729}}),
                "{value}"
            );
            let record_class = if flagged { json!("SNV") } else { Value::Null };
            assert_eq!(value["variant_class"], record_class, "{value}");

            // The text formats read the row's hash, which holds the record's
            // and the co-located value under the plugin's name once its flag
            // is on.
            let row = &v.transcript_consequences[1];
            let text =
                |field: &str| fields::field_value(field, &v, Some(row), &formatter.options, ",");
            let (frequency, old_frequency, significance, class) = if flagged {
                ("0.0004729", "0.25", "uncertain_significance", "SNV")
            } else {
                ("", "", "", "")
            };
            assert_eq!(text("gnomADg_AF"), frequency, "flags on: {flagged}");
            assert_eq!(text("ASN_AF"), old_frequency, "flags on: {flagged}");
            assert_eq!(text("CLIN_SIG"), significance, "flags on: {flagged}");
            assert_eq!(text("VARIANT_CLASS"), class, "flags on: {flagged}");
        }
    }

    /// A key of a `--custom` source that a plugin declares as a field is on an
    /// entry only as the plugin's own. The JSON factory moves the custom
    /// annotations to the record before any row hash is built
    /// (OutputFactory/JSON.pm:249-258), so no row's hash holds a source's
    /// names or its records' fields, where the text formats print them in the
    /// plugin's column (OutputFactory.pm:1300-1305).
    #[test]
    fn a_custom_key_a_plugin_declares_is_on_an_entry_only_where_the_plugin_wrote_it() {
        use crate::output::fields::{self, CustomSourceInfo};
        use vep_core::custom_annotation::{CustomAnnotations, CustomFieldValue, CustomRecord};
        let declared = ["SOURCE", "Comment", "Comment_FILTER"];
        let annotations = |short_name: &str| {
            CustomAnnotations::new(
                Arc::from(short_name),
                vec![CustomRecord::new(
                    "exact_1".to_string(),
                    None,
                    [(
                        "FILTER".to_string(),
                        CustomFieldValue::List(vec!["PASS".to_string()]),
                    )]
                    .into_iter()
                    .collect(),
                    None,
                )],
                true,
                None,
            )
        };
        let source = |short_name: &str| CustomSourceInfo {
            file: "/data/sites.vcf.gz".to_string(),
            short_name: short_name.to_string(),
            filter_field: true,
            ..Default::default()
        };
        let mut v = make_annotated_variant();
        v.custom_annotations = vec![annotations("SOURCE"), annotations("Comment")];
        let mut unwritten = v.transcript_consequences[0].clone();
        unwritten.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = "ENST00000307301".into();
            transcript
        }));
        for key in declared {
            v.transcript_consequences[0]
                .plugin_data
                .insert(key.to_string(), "from_the_table".to_string());
        }
        v.transcript_consequences.push(unwritten);

        let formatter = JsonOutputFormatter::new(
            FieldOptions {
                custom_sources: vec![source("SOURCE"), source("Comment")],
                ..Default::default()
            },
            declared.map(str::to_string).to_vec(),
            Some("GRCh38".into()),
        );
        let value = formatter.format_variant_value(&v);
        let record = json!({"fields": {"FILTER": ["PASS"]}, "name": "exact_1"});
        assert_eq!(
            value["custom_annotations"],
            json!({"SOURCE": [record.clone()], "Comment": [record]}),
            "{value}"
        );
        let entry_of = |transcript: &str| {
            value["transcript_consequences"]
                .as_array()
                .expect("the transcript entries")
                .iter()
                .find(|entry| entry["transcript_id"] == transcript)
                .unwrap_or_else(|| panic!("no entry of {transcript}: {value}"))
        };
        let written = entry_of("ENST00000352957");
        let entry = entry_of("ENST00000307301")
            .as_object()
            .expect("the entry the plugin left unwritten");
        for key in ["source", "comment", "comment_filter"] {
            assert_eq!(written[key], "from_the_table", "{value}");
            assert!(
                !entry.contains_key(key),
                "{key} on the entry the plugin left unwritten: {value}"
            );
        }

        // The text formats read the row's hash, which holds the source's
        // names and its records' `FILTER` under the plugin's names.
        let row = &v.transcript_consequences[1];
        let text = |field: &str| fields::field_value(field, &v, Some(row), &formatter.options, ",");
        assert_eq!(text("SOURCE"), "exact_1");
        assert_eq!(text("Comment"), "exact_1");
        assert_eq!(text("Comment_FILTER"), "PASS");
    }
}
