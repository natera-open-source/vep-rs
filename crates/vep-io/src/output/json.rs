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
    pub fn format_record_value(&self, record: &[&InputVariant]) -> Value {
        let mut root: IndexMap<String, Value> = IndexMap::new();
        let first = record[0];
        if record.iter().any(|v| v.vep_skip) {
            // VEP emits the unannotated record as its input line alone.
            if let Some(raw) = &first.raw_input {
                root.insert("input".to_string(), Value::String(raw.clone()));
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
            root.insert("input".to_string(), Value::String(raw.clone()));
        }
        if self.options.variant_class {
            root.insert(
                "variant_class".to_string(),
                Value::String(first.variant_class.so_term().to_string()),
            );
        }
        if !first.existing_variation.is_empty() {
            root.insert(
                "existing_variation".to_string(),
                Value::Array(
                    first
                        .existing_variation
                        .iter()
                        .map(|s| Value::String(s.clone()))
                        .collect(),
                ),
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
                .map(|cv| Value::Object(build_colocated_variant_json(cv)))
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
        let mut obj: IndexMap<String, Value> = IndexMap::new();
        let mut names: Vec<&str> = DEFAULT_OUTPUT_COLS
            .iter()
            .copied()
            .filter(|c| !matches!(*c, "Uploaded_variation" | "Location" | "Existing_variation"))
            .collect();
        for k in fields::extra_keys(variant, tc, &self.options) {
            if !names.contains(&k) {
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
                    let (s, e) = value.split_once('-').unwrap_or((&value, &value));
                    if let Ok(n) = s.parse::<u64>() {
                        obj.insert(format!("{coord}_start"), Value::Number(n.into()));
                    }
                    let e = if e.chars().all(|c| c.is_ascii_digit()) && !e.is_empty() {
                        e
                    } else {
                        s
                    };
                    if let Ok(n) = e.parse::<u64>() {
                        obj.insert(format!("{coord}_end"), Value::Number(n.into()));
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
                "FLAGS" | "MANE" | "CLIN_SIG" | "PUBMED" => {
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

/// Build a JSON object for a co-located variant.
fn build_colocated_variant_json(
    cv: &vep_core::variant::ColocatedVariant,
) -> serde_json::Map<String, Value> {
    let mut obj: IndexMap<String, Value> = IndexMap::new();
    obj.insert("id".to_string(), Value::String(cv.id.clone()));
    if let Some(ref as_str) = cv.allele_string {
        obj.insert("allele_string".to_string(), Value::String(as_str.clone()));
    }
    obj.insert("start".to_string(), Value::Number(cv.start.into()));
    obj.insert("end".to_string(), Value::Number(cv.end.into()));
    obj.insert(
        "strand".to_string(),
        Value::Number((cv.strand as i64).into()),
    );
    if !cv.frequencies.is_empty() {
        let mut freq_map = IndexMap::new();
        for (key, &val) in &cv.frequencies {
            if let Some(n) = serde_json::Number::from_f64(val) {
                freq_map.insert(key.clone(), Value::Number(n));
            }
        }
        if !freq_map.is_empty() {
            obj.insert(
                "frequencies".to_string(),
                Value::Object(serde_json::Map::from_iter(freq_map)),
            );
        }
    }
    if !cv.clin_sig.is_empty() {
        obj.insert(
            "clin_sig".to_string(),
            Value::Array(
                cv.clin_sig
                    .iter()
                    .map(|s| Value::String(s.clone()))
                    .collect(),
            ),
        );
    }
    if cv.somatic {
        obj.insert("somatic".to_string(), Value::Number(1.into()));
    }
    if cv.phenotype {
        obj.insert("phenotype_or_disease".to_string(), Value::Number(1.into()));
    }
    if !cv.pubmed.is_empty() {
        obj.insert(
            "pubmed".to_string(),
            Value::Array(cv.pubmed.iter().map(|s| Value::String(s.clone())).collect()),
        );
    }
    serde_json::Map::from_iter(obj)
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
        assert!(
            parsed.get("variant_class").is_none(),
            "variant_class needs its flag"
        );
        assert!(parsed.get("existing_variation").is_none());
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

    #[test]
    fn multi_allelic_record_is_one_object() {
        let mut a = make_annotated_variant();
        a.id = None;
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

    #[test]
    fn unknown_coordinate_is_dropped() {
        let mut v = make_annotated_variant();
        v.transcript_consequences[0].cds_position = Some("493-?".into());
        let parsed = fmt(FieldOptions::default()).format_variant_value(&v);
        let tc = &parsed["transcript_consequences"][0];
        assert_eq!(tc["cds_start"], 493);
        assert_eq!(tc["cds_end"], 493);
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
