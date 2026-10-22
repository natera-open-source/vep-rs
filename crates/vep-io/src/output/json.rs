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
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

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
    /// `assembly` is the name VEP writes as `assembly_name` (the `--assembly`
    /// flag, else the cache's); `None` omits the key.
    pub fn new(
        options: FieldOptions,
        plugin_fields: Vec<String>,
        assembly: Option<String>,
    ) -> Self {
        let vf_level_fields = fields::vf_level_flag_fields(&options);
        Self {
            options,
            vf_level_fields,
            plugin_fields,
            assembly,
        }
    }

    /// The JSON value for one input record (all alleles split from one line).
    /// Every object, an unannotated record's included, carries `vep_rs_version`.
    pub fn format_record_value(&self, record: &[&InputVariant]) -> Value {
        let mut root: IndexMap<String, Value> = IndexMap::new();
        root.insert(
            "vep_rs_version".to_string(),
            Value::String(fields::VEP_RS_VERSION.to_string()),
        );
        let first = record[0];
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
                _ => {
                    let value = fields::field_value(name, first, None, &self.options, ",");
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
                    all_terms.push(
                        variant
                            .most_severe_consequence
                            .unwrap_or(Consequence::IntergenicVariant),
                    );
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
        // row's hash (OutputFactory.pm:1702-1706) and into an intergenic row's
        // when `--shift_3prime 1` gave its insertion or deletion a shift length,
        // 0 included (:1991-1998), so those entries carry a `location`.
        let include_location = match tc {
            Some(_) => self.options.shift_genomic,
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
        // The co-located variant and frequency fields are stubbed out of the
        // consequence hash by the JSON factory (OutputFactory/JSON.pm:493-515) and
        // written under `colocated_variants` instead; `SV` is the record's key
        // (OutputFactory.pm:962-964), written at the top level as
        // `colocated_structural_variants` and in no consequence object.
        for k in fields::extra_keys(variant, tc, &self.options) {
            if !names.contains(&k)
                && !self.vf_level_fields.contains(&k)
                && !fields::colocated_field(k)
            {
                names.push(k);
            }
        }
        for p in &self.plugin_fields {
            if !names.contains(&p.as_str()) {
                names.push(p.as_str());
            }
        }
        for name in names {
            let value = fields::field_value(name, variant, tc, &self.options, ",");
            if value.is_empty() || (value == "-" && !fields::dash_is_a_value(name)) {
                continue;
            }
            match name {
                "Consequence" => {
                    obj.insert(
                        "consequence_terms".to_string(),
                        Value::Array(
                            value
                                .split(',')
                                .map(|t| Value::String(t.to_string()))
                                .collect(),
                        ),
                    );
                }
                "Feature_type" | "NEAREST" => {}
                "Feature" => {
                    let key = format!(
                        "{}_id",
                        tc.map(feature_type_key)
                            .unwrap_or_else(|| "intergenic".to_string())
                    );
                    obj.insert(key, Value::String(value));
                }
                "cDNA_position" | "CDS_position" | "Protein_position" => {
                    let coord = name.trim_end_matches("_position").to_lowercase();
                    let (start, end) = position_bounds(&value);
                    if let Some(start) = start {
                        obj.insert(format!("{coord}_start"), start);
                    }
                    if let Some(end) = end {
                        obj.insert(format!("{coord}_end"), end);
                    }
                }
                "SIFT" | "PolyPhen" => {
                    let tool = name.to_lowercase();
                    let (prediction, score) = parse_prediction_score(&value);
                    if let Some(prediction) = prediction {
                        obj.insert(
                            format!("{tool}_prediction"),
                            Value::String(prediction.to_string()),
                        );
                    }
                    if let Some(score) = score {
                        let key = format!("{tool}_score");
                        let number = numberify(&key, score.to_string());
                        obj.insert(key, number);
                    }
                }
                "DOMAINS" => {
                    let domains =
                        domains_json(tc.map(|x| x.domains.as_slice()).unwrap_or_default());
                    obj.insert("domains".to_string(), Value::Array(domains));
                }
                // The fields VEP's output hash holds as lists (OutputFactory.pm:1478,
                // 1540, 1572, 1592-1600, 1662 at release/116.2), encoded as arrays;
                // the one undefined element VEP can produce, the miRNA structure
                // read one past its end, is `null`, while an empty element of a
                // split accession list stays an empty string.
                "FLAGS" | "MANE" | "RefSeq" | "SWISSPROT" | "TREMBL" | "UNIPARC"
                | "UNIPROT_ISOFORM" | "miRNA" => {
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
                _ => {
                    let key = json_key(name);
                    let v = if value == "YES" {
                        Value::Number(1.into())
                    } else {
                        numberify(&key, value)
                    };
                    obj.insert(key, v);
                }
            }
        }
        serde_json::Map::from_iter(obj)
    }
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
/// `--total_length` suffix (`N/T`) stays a string and is written at both ends;
/// the end of `N-M/T` falls back to `N`, `N-?/T` writes the start alone and
/// `?-N/T` writes neither.
fn position_bounds(value: &str) -> (Option<Value>, Option<Value>) {
    let (start, end) = value.split_once('-').unwrap_or((value, value));
    let end = if end.contains('?') || is_digits(end) {
        end
    } else {
        start
    };
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

/// VEP's JSON key for a field name: lowercased, then renamed where VEP renames.
fn json_key(field: &str) -> String {
    match field {
        "Gene" => "gene_id".to_string(),
        "Allele" => "variant_allele".to_string(),
        "SYMBOL" => "gene_symbol".to_string(),
        "SYMBOL_SOURCE" => "gene_symbol_source".to_string(),
        "OverlapBP" => "bp_overlap".to_string(),
        "OverlapPC" => "percentage_overlap".to_string(),
        "RefSeq" => "refseq_transcript_ids".to_string(),
        "ENSP" => "protein_id".to_string(),
        "SV" => "colocated_structural_variants".to_string(),
        other => other.to_lowercase(),
    }
}

/// Numeric-looking values become numbers, as VEP's `numberify` does, except the
/// identifier keys it exempts.
fn numberify(key: &str, value: String) -> Value {
    if matches!(
        key,
        "seq_region_name" | "id" | "gene_id" | "gene_symbol" | "transcript_id"
    ) {
        return Value::String(value);
    }
    if let Ok(i) = value.parse::<i64>() {
        return Value::Number(i.into());
    }
    if value.contains('.') || value.contains('e') || value.contains('E') {
        if let Ok(f) = value.parse::<f64>() {
            // Perl's numeric conversion drops a zero fraction: "100.00" is 100.
            if f.fract() == 0.0 && f.abs() < 1e15 {
                return Value::Number((f as i64).into());
            }
            if let Some(n) = serde_json::Number::from_f64(f) {
                return Value::Number(n);
            }
        }
    }
    Value::String(value)
}

/// Splits a `SIFT` or `PolyPhen` value into its prediction and score as
/// `OutputFactory/JSON.pm:353-360` (release/116.2) does with
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
/// (release/116.2) builds it: each entry the `DOMAINS` column prints (the label
/// and identifier parts `grep {$_}` keeps, joined by `:`, OutputFactory.pm:1500)
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
/// (OutputFactory/JSON.pm:407-476 at release/116.2) writes the record's hash:
/// every cache column that is set, not empty and not `0` (:433-435), the
/// `failed` and matched-allele bookkeeping dropped, `chr` as `seq_region_name`
/// and `variation_name` as `id`
/// (:98-112), `clin_sig` and `pubmed` as lists, `var_synonyms` as
/// `{source: [ids]}` split on `--`, `::` and `,` (:449-471), `clinical_impact` as
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
    if !cv.clin_sig.is_empty() {
        obj.insert("clin_sig".to_string(), numberified_list(&cv.clin_sig));
    }
    if !cv.pubmed.is_empty() {
        obj.insert("pubmed".to_string(), numberified_list(&cv.pubmed));
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
    if let Some(synonyms) = cv.var_synonyms.as_deref().filter(|s| perl_true(s)) {
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
/// exemption).
fn numberified_list(values: &[String]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|value| numberify("", value.clone()))
            .collect(),
    )
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
    use serde_json::json;
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, FeatureType, Impact, TranscriptConsequence};
    use vep_core::terms::Vocabulary;

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
        let tc = TranscriptConsequence {
            transcript_id: "ENST00000352957".into(),
            gene_id: "ENSG00000141959".into(),
            gene_symbol: Some("MRPL39".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:14027".into()),
            consequences: smallvec![Consequence::SynonymousVariant],
            impact: Impact::LOW,
            biotype: Some("protein_coding".into()),
            canonical: true,
            cdna_position: Some("428".into()),
            cds_position: Some("366".into()),
            protein_position: Some("122".into()),
            amino_acids: Some("L".into()),
            codons: Some("ctG/ctA".into()),
            strand: -1,
            feature_type: FeatureType::Transcript,
            flags: std::sync::Arc::from(["cds_end_NF".to_string()]),
            ..Default::default()
        };
        variant.transcript_consequences = vec![tc];
        variant.most_severe_consequence = Some(Consequence::SynonymousVariant);
        variant
    }

    fn fmt(options: FieldOptions) -> JsonOutputFormatter {
        JsonOutputFormatter::new(options, Vec::new(), Some("GRCh37".into()))
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
    /// Ensembl VEP 116.2 writes with a database for `21 10538885 A C`
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
        v.colocated_variants = vec![vep_core::variant::ColocatedVariant {
            id: "rs123".into(),
            ..Default::default()
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
    /// (release 116.2), which keeps a `?` end as the end and then deletes any end
    /// that is not a number, and the row of the same shape Ensembl VEP 116.2
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
    /// follows ensembl-vep `OutputFactory/JSON.pm:337-341` (release 116.2), which
    /// deletes any start that is not a number, and the row of the same shape
    /// Ensembl VEP 116.2 writes for `tests/golden/116/GRCh38` record `4287535`
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

    /// `OutputFactory/JSON.pm` on the `--total_length` forms: `N/T` is kept as a
    /// string at both ends (`numberify` leaves it alone); the end of `N-M/T` falls
    /// back to the start `N`; `N-?/T` writes the start alone; `?-N/T` writes
    /// neither.
    #[test]
    fn total_length_forms_follow_the_perl_split() {
        assert_eq!(
            position_bounds("428/1234"),
            (Some(json!("428/1234")), Some(json!("428/1234")))
        );
        assert_eq!(
            position_bounds("428-430/1234"),
            (Some(json!(428)), Some(json!(428)))
        );
        assert_eq!(position_bounds("428-?/1234"), (Some(json!(428)), None));
        assert_eq!(position_bounds("?-430/1234"), (None, None));
    }

    /// The `--total_length` columns as the JSON writer files them, Ensembl VEP's
    /// shape (`OutputFactory/JSON.pm` `add_VariationFeatureOverlapAllele_info`): a
    /// single position with its total is the string `"53/500"` at both ends, a
    /// range with its total keeps its start as a number at both ends because the
    /// `/total` tail is not a plain number, and `consequence_terms` is untouched.
    #[test]
    fn total_length_json_keeps_perls_split_strings() {
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
        assert_eq!((&tc["cds_start"], &tc["cds_end"]), (&json!(4), &json!(4)));
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
    /// `existing_variation`. Ensembl VEP 116.2 writes
    /// `{"var_synonyms": {"ClinVar": ["RCV000886041", "VCV000713837"]}, "clin_sig":
    /// ["likely_benign"], "clin_sig_allele": "T:likely_benign", "clin_sig_ref_allele":
    /// "G", "phenotype_or_disease": 1, "id": "rs374325784", "seq_region_name": "21",
    /// "start": 10541121, "end": 10541121, "strand": 1, "allele_string": "G/A/C/T"}`
    /// for `21 10541121 G A` of the `GRCh38-colocated-alt` corpus.
    #[test]
    fn var_synonyms_nested_by_source() {
        let mut v = make_annotated_variant();
        v.colocated_variants = vec![vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs374325784".into(),
            start: 10_541_121,
            end: 10_541_121,
            allele_string: Some("G/A/C/T".into()),
            strand: 1,
            clin_sig: vec!["likely_benign".into()],
            phenotype: true,
            pubmed: vec!["26076356".into(), "11524731".into()],
            clinvar: Some(Box::new(vep_core::variation::ClinVarDetail {
                clin_sig_allele: Some("T:likely_benign".into()),
                clin_sig_ref_allele: Some("G".into()),
                clinical_impact: Some("Neoplasm(Tier I - Strong)".into()),
            })),
            var_synonyms: Some(
                "ClinVar::RCV000886041,VCV000713837--OMIM::605717.0001--UniProt::VAR_089699".into(),
            ),
            ..Default::default()
        }];
        v.existing_variation = vec!["rs374325784".into()];
        let options = FieldOptions {
            check_existing: true,
            pubmed: true,
            ..Default::default()
        };
        let parsed = fmt(options).format_variant_value(&v);
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
        // A germline record without ClinVar fields or synonyms carries the
        // coordinates, alleles and name alone (:433-435 removes the rest).
        v.colocated_variants = vec![vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs1300916113".into(),
            start: 5_030_602,
            end: 5_030_602,
            allele_string: Some("A/C".into()),
            strand: 1,
            ..Default::default()
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
    /// Ensembl VEP 116.2 writes `minor_allele` and `minor_allele_freq` in every
    /// object on a cache declaring the two pre-108 columns; a release-116 cache
    /// declares no such column.
    #[test]
    fn colocated_variant_echoes_the_other_cache_columns() {
        let mut v = make_annotated_variant();
        v.colocated_variants = vec![vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs1".into(),
            start: 1_004,
            end: 1_004,
            allele_string: Some("C/A".into()),
            strand: 1,
            frequency_strings: [
                ("minor_allele", "A"),
                ("minor_allele_freq", "0.12"),
                ("AF", "A:0.5"),
                ("other_af", "0"),
                ("some_count", "0"),
                ("some_text", "x"),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
            ..Default::default()
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
    /// (:160-163): Ensembl VEP 116.2 writes
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
        let rs1 = vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs1".into(),
            start: 25000100,
            end: 25000100,
            allele_string: Some("A/G/T".into()),
            strand: 1,
            matched_alleles: vec![matched("G", 0), matched("T", 1)],
            frequency_strings: [
                ("AF", "G:1"),
                ("gnomADg_AFR", "T:0,G:1"),
                ("gnomADg", "T:0.0002047,G:0.9994"),
            ]
            .into_iter()
            .map(|(column, value)| (column.to_string(), value.to_string()))
            .collect(),
            ..Default::default()
        };
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
    /// written spelling and finds no entry: Ensembl VEP 116.2 writes
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
        let rs = vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs1370167485".into(),
            start: 5030602,
            end: 5030602,
            allele_string: Some("A/AA".into()),
            strand: 1,
            matched_alleles: vec![matched],
            frequency_strings: [("gnomADg".to_string(), "AA:0.0001735".to_string())].into(),
            ..Default::default()
        };
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
    /// Ensembl VEP's JSON (`[split(',', $id)]`, OutputFactory.pm:1572 at
    /// release/116.2, encoded as is by `JSON.pm`), one element per
    /// comma-separated accession, and a `-` value sets no key: the GRCh37 APP
    /// transcript's `TREMBL=L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN` is a
    /// three-element array.
    #[test]
    fn uniprot_columns_are_arrays() {
        let tc = attribute_object(|tc| {
            tc.swissprot = Some("A4_HUMAN".to_string());
            tc.trembl = Some("L7XE61_HUMAN,L7XCZ9_HUMAN,B7Z313_HUMAN".to_string());
            tc.uniparc = Some("UPI000002DB1C".to_string());
            tc.uniprot_isoform = Some("-".to_string());
        });
        assert_eq!(tc["swissprot"], json!(["A4_HUMAN"]));
        assert_eq!(
            tc["trembl"],
            json!(["L7XE61_HUMAN", "L7XCZ9_HUMAN", "B7Z313_HUMAN"])
        );
        assert_eq!(tc["uniparc"], json!(["UPI000002DB1C"]));
        assert!(tc.get("uniprot_isoform").is_none(), "{tc}");
        let empty_element = attribute_object(|tc| tc.trembl = Some("A,,B".to_string()));
        assert_eq!(empty_element["trembl"], json!(["A", "", "B"]));
    }

    /// `RefSeq` is renamed `refseq_transcript_ids` (JSON.pm:106) and is the
    /// array of the comma-separated ids (OutputFactory.pm:1540):
    /// `RefSeq=NM_198996.4,NM_001302998.2` on the LIPI transcript
    /// ENST00000681601 is a two-element array.
    #[test]
    fn refseq_transcript_ids_is_array() {
        let tc = attribute_object(|tc| {
            tc.refseq = Some("NM_198996.4,NM_001302998.2".to_string());
        });
        assert_eq!(
            tc["refseq_transcript_ids"],
            json!(["NM_198996.4", "NM_001302998.2"])
        );
        assert!(tc.get("refseq").is_none());
        let dash = attribute_object(|tc| tc.refseq = Some("-".to_string()));
        assert!(dash.get("refseq_transcript_ids").is_none());
    }

    /// `mirna` is the array of structure terms (OutputFactory.pm:1662), `mane`
    /// the array of set names (:1592, :1600); an undefined element is `null`.
    #[test]
    fn mirna_and_mane_are_arrays() {
        let tc = attribute_object(|tc| {
            tc.mirna = vec!["miRNA_loop".to_string(), "miRNA_stem".to_string()];
            tc.mane_select = Some("NM_021075.4".to_string());
            tc.mane_plus_clinical = Some("NM_058174.3".to_string());
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
    /// 116.2 writes `["nsv2770691","nsv3109934","nsv3163927","nsv3167562",
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
    /// `$1` and `$2` keeping the last successful match: the Ensembl VEP 116.2
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
}
