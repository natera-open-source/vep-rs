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

use super::fields::{self, FieldOptions, DEFAULT_OUTPUT_COLS};
use super::OutputFormatter;

/// JSON output formatter.
pub struct JsonOutputFormatter {
    options: FieldOptions,
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
        Self {
            options,
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
        root.insert(
            "allele_string".to_string(),
            Value::String(string_source.record_allele_string()),
        );
        if let Some(asm) = &self.assembly {
            root.insert("assembly_name".to_string(), Value::String(asm.clone()));
        }
        if let Some(raw) = &first.raw_input {
            root.insert("input".to_string(), Value::String(raw.to_string()));
        }
        if self.options.variant_class {
            root.insert(
                "variant_class".to_string(),
                Value::String(first.variant_class.so_term().to_string()),
            );
        }
        let mut buckets: IndexMap<String, Vec<Value>> = IndexMap::new();
        let mut all_terms: Vec<Consequence> = Vec::new();
        for (variant, tc) in fields::record_rows(record) {
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
            .map(|c| c.so_term().to_string())
            .unwrap_or_else(|| "?".to_string());
        root.insert(
            "most_severe_consequence".to_string(),
            Value::String(most_severe),
        );

        if !first.colocated_variants.is_empty() {
            let cvs: Vec<Value> = first
                .colocated_variants
                .iter()
                .map(|cv| Value::Object(colocated_variant_object(cv, record)))
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
            let consequence_terms = terms
                .iter()
                .map(|c| Value::String(c.so_term().to_string()))
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
        // written under `colocated_variants` instead.
        for k in fields::extra_keys(variant, tc, &self.options) {
            if !names.contains(&k) && !fields::colocated_field(k) {
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
            if value.is_empty() || (value == "-" && name != "Allele") {
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
                "Feature_type" => {}
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
                    let (pred, score) = parse_prediction_score(&value);
                    if let Some(p) = pred.filter(|p| !p.is_empty()) {
                        obj.insert(format!("{tool}_prediction"), Value::String(p));
                    }
                    if let Some(s) = score.and_then(serde_json::Number::from_f64) {
                        obj.insert(format!("{tool}_score"), Value::Number(s));
                    }
                }
                "DOMAINS" => {
                    let doms: Vec<Value> = value
                        .split(',')
                        .filter_map(|d| d.split_once(':'))
                        .map(|(db, id)| {
                            Value::Object(serde_json::Map::from_iter([
                                ("db".to_string(), Value::String(db.to_string())),
                                ("name".to_string(), Value::String(id.to_string())),
                            ]))
                        })
                        .collect();
                    obj.insert("domains".to_string(), Value::Array(doms));
                }
                "FLAGS" | "MANE" => {
                    obj.insert(
                        json_key(name),
                        Value::Array(
                            value
                                .split(',')
                                .map(|t| Value::String(t.to_string()))
                                .collect(),
                        ),
                    );
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

/// Parse a SIFT/PolyPhen string like "deleterious(0.01)" into (prediction, score).
fn parse_prediction_score(value: &str) -> (Option<String>, Option<f64>) {
    if let Some(open_paren) = value.find('(') {
        if let Some(close_paren) = value.find(')') {
            let prediction = value[..open_paren].to_string();
            let score = value[open_paren + 1..close_paren].parse::<f64>().ok();
            return (Some(prediction), score);
        }
    }
    (Some(value.to_string()), value.parse::<f64>().ok())
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
/// the raw string under `clinvar_somatic_classification`, and the frequency
/// columns replaced by `frequencies`: for each alternate allele of the record
/// the values of the record allele matched to it, keyed by the lower-cased
/// population (:419-431). Every value is numberified except the identifiers
/// (:288, `%NUMBERIFY_EXEMPT`).
///
/// `record` is every allele split from the input line; the record's copy of the
/// variant on allele `a` carries the frequencies matched to `a`.
fn colocated_variant_object(
    cv: &vep_core::variant::ColocatedVariant,
    record: &[&InputVariant],
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
    let mut frequencies: IndexMap<String, Value> = IndexMap::new();
    for allele in record {
        let copy = allele
            .colocated_variants
            .iter()
            .find(|c| c.id == cv.id && c.start == cv.start && c.end == cv.end);
        let Some(copy) = copy.filter(|c| !c.frequencies.is_empty()) else {
            continue;
        };
        let mut by_population: IndexMap<String, Value> = IndexMap::new();
        for population in FREQ_FIELDS {
            if let Some(&value) = copy.frequencies.get(*population) {
                by_population.insert(
                    population.to_lowercase(),
                    numberify("frequency", value.to_string()),
                );
            }
        }
        if !by_population.is_empty() {
            frequencies.insert(
                allele.display_allele(),
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

    /// `OutputFactory/JSON.pm:333-343` (release 116.2) on `--total_length` values:
    /// a point `428/1200` is split on `-` into one part, so `cdna_start` and
    /// `cdna_end` are both the string `"428/1200"` (`numberify` leaves a value with
    /// a `/` alone); a range `900-901/1200` has an end that fails `^\d+$`, so
    /// `cds_end` falls back to the start `900` while `cds_start` is the number.
    #[test]
    fn total_length_json_keeps_perls_split_strings() {
        let tc = positions("428/1200", "900-901/1200", "122/400");
        assert_eq!(tc["cdna_start"], json!("428/1200"));
        assert_eq!(tc["cdna_end"], json!("428/1200"));
        assert_eq!(tc["cds_start"], json!(900));
        assert_eq!(tc["cds_end"], json!(900));
        assert_eq!(tc["protein_start"], json!("122/400"));
        assert_eq!(tc["protein_end"], json!("122/400"));
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

    /// `frequencies` nests, under each alternate allele of the record, the
    /// matched record allele's values keyed by the lower-cased population
    /// (OutputFactory/JSON.pm:419-431): an allele the variant was not matched to
    /// has no entry, and a record without matched frequencies carries no
    /// `frequencies` key.
    #[test]
    fn colocated_frequencies_nest_per_allele() {
        let mut a = make_annotated_variant();
        a.allele_index = 0;
        a.record_allele_string_multi = Some("A/G/T".into());
        let mut g_copy = vep_core::variant::ColocatedVariant {
            chr: Some("21".into()),
            id: "rs1".into(),
            start: 25000100,
            end: 25000100,
            allele_string: Some("A/G/T".into()),
            strand: 1,
            ..Default::default()
        };
        g_copy
            .frequencies
            .insert("gnomADg_AFR".to_string(), 6.875e-5);
        g_copy.frequencies.insert("AF".to_string(), 1.0);
        a.colocated_variants = vec![g_copy.clone()];
        let mut b = a.clone();
        b.alt_alleles = vec![b"T".to_vec()];
        b.allele_string = "A/T".into();
        b.allele_index = 1;
        let mut t_copy = g_copy.clone();
        t_copy.frequencies.clear();
        t_copy.frequencies.insert("gnomADg_AFR".to_string(), 0.0);
        b.colocated_variants = vec![t_copy];
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert_eq!(
            parsed["colocated_variants"][0]["frequencies"],
            json!({
                "G": {"af": 1, "gnomadg_afr": 6.875e-5},
                "T": {"gnomadg_afr": 0}
            })
        );
        b.colocated_variants[0].frequencies.clear();
        a.colocated_variants[0].frequencies.clear();
        let parsed = fmt(FieldOptions::default()).format_record_value(&[&a, &b]);
        assert!(parsed["colocated_variants"][0].get("frequencies").is_none());
    }

    #[test]
    fn no_header_for_json() {
        let mut f = fmt(FieldOptions::default());
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn prediction_score_parsing() {
        assert_eq!(
            parse_prediction_score("deleterious(0.01)"),
            (Some("deleterious".into()), Some(0.01))
        );
        assert_eq!(
            parse_prediction_score("tolerated"),
            (Some("tolerated".into()), None)
        );
    }
}
