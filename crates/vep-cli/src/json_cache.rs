// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Load transcripts from the converted JSON cache directory.
//!
//! The JSON cache (written by `vep-cache-builder` or converted from a Perl VEP cache) uses field names that
//! differ from the Rust `Transcript` struct (e.g. `dbID` vs `db_id`,
//! `is_canonical` vs `canonical`). This module provides deserialization types
//! that map the JSON format and converts them to the canonical Rust types.
//!
//! Perl citations name modules of ensembl-vep release/115 (`Bio/EnsEMBL/VEP/...`);
//! `VariationEffect.pm` is `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm` in
//! ensembl-variation release/115.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Deserialize;
use tracing::info;

use vep_core::coordinate::Strand;
use vep_core::transcript::{
    Attribute, Exon, ExonCoordMapper, Intron, MapperPair, PredictionMatrix, ProteinFeature,
    ProteinFunctionPredictions, Transcript, TranscriptMapper, TranscriptVEFC, Translation,
};
use vep_core::variation::CachedVariation;

/// Visitors for the numeric fields a Storable-derived cache may store as a JSON number
/// or as a string holding one.
///
/// An `#[serde(untagged)]` enum would buffer each value and then try its
/// variants in turn, formatting an error message for every variant that does
/// not match, so a plain integer field cost an allocation and a formatted
/// error on the load path. Each visitor below accepts the same inputs with the
/// same results: an integer as itself, a negative integer as 0, a float
/// truncated, a string parsed (an empty string is "absent" where the field is
/// optional), and for `u8` a bool as 0/1.
macro_rules! num_or_str_visitor {
    ($name:ident, $t:ty) => {
        num_or_str_visitor!(@impl $name, $t,);
    };
    ($name:ident, $t:ty, bool) => {
        num_or_str_visitor!(
            @impl $name, $t,
            fn visit_bool<E: serde::de::Error>(self, v: bool) -> Result<Self::Value, E> {
                Ok(Some(v as $t))
            }
        );
    };
    (@impl $name:ident, $t:ty, $($extra:tt)*) => {
        struct $name;
        impl serde::de::Visitor<'_> for $name {
            type Value = Option<$t>;

            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str(concat!(
                    "a ",
                    stringify!($t),
                    ", a string holding one, or an empty string"
                ))
            }

            fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
                Ok(Some(v as $t))
            }

            fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
                Ok(Some(if v < 0 { 0 } else { v as $t }))
            }

            fn visit_f64<E: serde::de::Error>(self, v: f64) -> Result<Self::Value, E> {
                Ok(Some(v as $t))
            }

            $($extra)*

            fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
                if v.is_empty() {
                    Ok(None)
                } else {
                    v.parse::<$t>().map(Some).map_err(E::custom)
                }
            }
        }
    };
}

num_or_str_visitor!(U64OrStr, u64);
num_or_str_visitor!(U32OrStr, u32);
num_or_str_visitor!(U8OrStr, u8, bool);

/// `Option<T>` through `deserialize_option`, so a JSON `null` is `None` and
/// anything else goes to `inner`.
struct OptionOf<V>(V);
impl<'de, T, V> serde::de::Visitor<'de> for OptionOf<V>
where
    V: serde::de::Visitor<'de, Value = Option<T>>,
{
    type Value = Option<T>;

    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0.expecting(f)?;
        f.write_str(", or null")
    }

    fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }

    fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
        Ok(None)
    }

    fn visit_some<D: serde::Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        d.deserialize_any(self.0)
    }
}

/// Deserialize helper: parse a value that may be a JSON number or a string
/// containing a number.
fn deserialize_string_or_number<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer
        .deserialize_any(U64OrStr)?
        .ok_or_else(|| serde::de::Error::custom("cannot parse integer from empty string"))
}

fn deserialize_optional_string_or_number<'de, D>(deserializer: D) -> Result<Option<u64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer.deserialize_option(OptionOf(U64OrStr))
}

/// Visitor for a signed field stored as a JSON integer or a string; a string
/// that does not parse yields `fallback`.
struct I64OrStr {
    fallback: i64,
}
impl serde::de::Visitor<'_> for I64OrStr {
    type Value = i64;

    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("an integer or a string")
    }

    fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
        Ok(v)
    }

    fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
        i64::try_from(v).map_err(|_| E::invalid_value(serde::de::Unexpected::Unsigned(v), &self))
    }

    fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
        Ok(v.parse::<i64>().unwrap_or(self.fallback))
    }
}

fn deserialize_i8_from_number<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer
        .deserialize_any(I64OrStr { fallback: 0 })
        .map(|n| n as i8)
}

fn deserialize_optional_u8<'de, D>(deserializer: D) -> Result<Option<u8>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer.deserialize_option(OptionOf(U8OrStr))
}

fn deserialize_optional_u32<'de, D>(deserializer: D) -> Result<Option<u32>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer.deserialize_option(OptionOf(U32OrStr))
}

/// Deserialize a usize that may be null in Storable-derived caches.
fn deserialize_optional_usize_as_zero<'de, D>(deserializer: D) -> Result<usize, D::Error>
where
    D: serde::Deserializer<'de>,
{
    Option::<usize>::deserialize(deserializer).map(|v| v.unwrap_or(0))
}

/// Deserialize a string straight into the shared `Arc<str>` the `Transcript`
/// field holds, without an intermediate `String`.
fn deserialize_arc_str<'de, D>(deserializer: D) -> Result<Arc<str>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct ArcStrVisitor;
    impl serde::de::Visitor<'_> for ArcStrVisitor {
        type Value = Arc<str>;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a string")
        }

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
            Ok(Arc::from(v))
        }

        fn visit_string<E: serde::de::Error>(self, v: String) -> Result<Self::Value, E> {
            Ok(Arc::from(v))
        }
    }
    deserializer.deserialize_str(ArcStrVisitor)
}

fn deserialize_optional_arc_str<'de, D>(deserializer: D) -> Result<Option<Arc<str>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct OptVisitor;
    impl<'de> serde::de::Visitor<'de> for OptVisitor {
        type Value = Option<Arc<str>>;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a string or null")
        }

        fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_some<Inner: serde::Deserializer<'de>>(
            self,
            deserializer: Inner,
        ) -> Result<Self::Value, Inner::Error> {
            deserialize_arc_str(deserializer).map(Some)
        }
    }
    deserializer.deserialize_option(OptVisitor)
}

fn deserialize_bool_from_int<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct BoolOrInt;
    impl serde::de::Visitor<'_> for BoolOrInt {
        type Value = bool;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a bool or an integer")
        }

        fn visit_bool<E: serde::de::Error>(self, v: bool) -> Result<Self::Value, E> {
            Ok(v)
        }

        fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
            Ok(v != 0)
        }

        fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
            i64::try_from(v)
                .map(|i| i != 0)
                .map_err(|_| E::invalid_value(serde::de::Unexpected::Unsigned(v), &self))
        }
    }
    deserializer.deserialize_any(BoolOrInt)
}

#[derive(Deserialize)]
struct JsonTranscript {
    #[serde(deserialize_with = "deserialize_arc_str")]
    stable_id: Arc<str>,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    version: Option<u32>,
    #[serde(
        rename = "dbID",
        default,
        deserialize_with = "deserialize_optional_string_or_number"
    )]
    db_id: Option<u64>,
    #[serde(deserialize_with = "deserialize_arc_str")]
    gene_stable_id: Arc<str>,
    #[serde(default)]
    chr: Option<String>,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(deserialize_with = "deserialize_strand")]
    strand: Strand,
    #[serde(deserialize_with = "deserialize_arc_str")]
    biotype: Arc<str>,
    source: String,
    #[serde(default)]
    description: Option<String>,
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_symbol: Option<Arc<str>>,
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_symbol_source: Option<Arc<str>>,
    #[serde(
        default,
        rename = "gene_hgnc_id",
        deserialize_with = "deserialize_optional_arc_str"
    )]
    hgnc_id: Option<Arc<str>>,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    gene_phenotype: Option<u32>,
    #[serde(
        default,
        rename = "is_canonical",
        deserialize_with = "deserialize_bool_from_int"
    )]
    canonical: bool,
    #[serde(default)]
    mane_select: Option<String>,
    #[serde(default)]
    mane_plus_clinical: Option<String>,
    #[serde(default)]
    tsl: Option<u8>,
    #[serde(default)]
    appris: Option<String>,
    #[serde(default)]
    ccds: Option<String>,
    #[serde(
        default,
        rename = "protein",
        deserialize_with = "deserialize_optional_arc_str"
    )]
    protein_id: Option<Arc<str>>,
    #[serde(default)]
    refseq: Option<String>,
    #[serde(default)]
    swissprot: Option<String>,
    #[serde(default)]
    trembl: Option<String>,
    #[serde(default)]
    uniparc: Option<String>,
    #[serde(default)]
    exons: Vec<JsonExon>,
    #[serde(default, deserialize_with = "deserialize_optional_string_or_number")]
    cdna_coding_start: Option<u64>,
    #[serde(default, deserialize_with = "deserialize_optional_string_or_number")]
    cdna_coding_end: Option<u64>,
    #[serde(default)]
    translation: Option<JsonTranslation>,
    #[serde(default)]
    attributes: Vec<JsonAttribute>,
    #[serde(default, rename = "variation_effect_feature_cache")]
    vefc: Option<JsonVEFC>,
}

fn deserialize_strand<'de, D>(deserializer: D) -> Result<Strand, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let n = deserializer.deserialize_any(I64OrStr { fallback: 1 })?;
    Ok(if n >= 0 {
        Strand::Forward
    } else {
        Strand::Reverse
    })
}

#[derive(Deserialize)]
struct JsonExon {
    #[serde(default)]
    stable_id: Option<String>,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    rank: Option<u32>,
    #[serde(default, deserialize_with = "deserialize_i8_from_number")]
    phase: i8,
    #[serde(default, deserialize_with = "deserialize_i8_from_number")]
    end_phase: i8,
}

#[derive(Deserialize)]
struct JsonIntron {
    #[serde(deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    rank: Option<u32>,
}

#[derive(Deserialize)]
struct JsonTranslation {
    stable_id: String,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    version: Option<u32>,
    #[serde(
        rename = "dbID",
        default,
        deserialize_with = "deserialize_optional_string_or_number"
    )]
    db_id: Option<u64>,
    #[serde(default, deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(default, deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(default)]
    seq: Option<String>,
}

#[derive(Deserialize)]
struct JsonAttribute {
    code: String,
    value: String,
}

#[derive(Deserialize)]
struct JsonVEFC {
    #[serde(default, deserialize_with = "deserialize_optional_u8")]
    codon_table: Option<u8>,
    #[serde(default)]
    five_prime_utr: Option<String>,
    #[serde(default)]
    three_prime_utr: Option<String>,
    #[serde(default)]
    translateable_seq: Option<String>,
    #[serde(default)]
    peptide: Option<String>,
    #[serde(default)]
    introns: Vec<JsonIntron>,
    #[serde(default)]
    sorted_exons: Vec<JsonExon>,
    #[serde(default)]
    mapper: Option<JsonMapper>,
    #[serde(default)]
    protein_features: Vec<JsonProteinFeature>,
    #[serde(default)]
    protein_function_predictions: Option<JsonProteinFunctionPredictions>,
}

/// JSON-serialized SIFT/PolyPhen prediction matrices.
///
/// Each matrix is stored as a base64-encoded gzipped binary blob.
/// See [`vep_core::prediction`] for format documentation.
#[derive(Deserialize)]
struct JsonProteinFunctionPredictions {
    #[serde(default)]
    sift: Option<JsonPredictionMatrix>,
    #[serde(default)]
    polyphen_humdiv: Option<JsonPredictionMatrix>,
    #[serde(default)]
    polyphen_humvar: Option<JsonPredictionMatrix>,
}

#[derive(Deserialize)]
struct JsonPredictionMatrix {
    #[serde(default)]
    analysis: String,
    #[serde(default)]
    sub_analysis: Option<String>,
    #[serde(default, deserialize_with = "deserialize_optional_usize_as_zero")]
    peptide_length: usize,
    /// Base64-encoded gzipped prediction matrix.
    #[serde(default)]
    predictions_data: String,
}

#[derive(Deserialize)]
struct JsonMapper {
    #[serde(default, deserialize_with = "deserialize_i8_from_number")]
    start_phase: i8,
    #[serde(default, deserialize_with = "deserialize_optional_string_or_number")]
    cdna_coding_start: Option<u64>,
    #[serde(default, deserialize_with = "deserialize_optional_string_or_number")]
    cdna_coding_end: Option<u64>,
    // `pair_count` is intentionally not deserialized: `ExonCoordMapper::new`
    // derives the count from `pairs.len()`. Serde ignores the JSON key.
    #[serde(default)]
    pairs: Vec<JsonMapperPair>,
}

#[derive(Deserialize)]
struct JsonMapperPair {
    #[serde(deserialize_with = "deserialize_string_or_number")]
    from_start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    from_end: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    to_start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    to_end: u64,
    #[serde(deserialize_with = "deserialize_i8_from_number")]
    ori: i8,
}

#[derive(Deserialize)]
struct JsonProteinFeature {
    #[serde(deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    hseqname: String,
    #[serde(default, deserialize_with = "deserialize_analysis")]
    analysis: Option<String>,
}

/// Deserialize `analysis` which may be a string (vep-cache-builder / GRCh37 Storable)
/// or an object with `display_label` (raw Storable conversion).
fn deserialize_analysis<'de, D>(deserializer: D) -> Result<Option<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de;

    struct AnalysisVisitor;
    impl<'de> de::Visitor<'de> for AnalysisVisitor {
        type Value = Option<String>;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a string, null, or an object with display_label")
        }

        fn visit_none<E: de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }

        fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
            if v.is_empty() {
                Ok(None)
            } else {
                Ok(Some(v.to_string()))
            }
        }

        fn visit_string<E: de::Error>(self, v: String) -> Result<Self::Value, E> {
            if v.is_empty() {
                Ok(None)
            } else {
                Ok(Some(v))
            }
        }

        fn visit_map<A: de::MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
            let mut display_label: Option<String> = None;
            while let Some(key) = map.next_key::<String>()? {
                let val: serde_json::Value = map.next_value()?;
                if key == "display_label" {
                    display_label = val.as_str().map(|s| s.to_string());
                }
            }
            Ok(display_label)
        }
    }

    deserializer.deserialize_any(AnalysisVisitor)
}

/// Decode a base64-encoded gzipped prediction matrix into raw bytes.
fn decode_prediction_matrix(jpm: JsonPredictionMatrix) -> Option<PredictionMatrix> {
    use base64::Engine;

    let compressed = base64::engine::general_purpose::STANDARD
        .decode(&jpm.predictions_data)
        .ok()?;
    let data = vep_core::prediction::gunzip_matrix(&compressed).ok()?;
    Some(PredictionMatrix {
        analysis: jpm.analysis,
        sub_analysis: jpm.sub_analysis,
        peptide_length: jpm.peptide_length,
        predictions_data: data,
    })
}

fn convert_protein_function_predictions(
    jp: JsonProteinFunctionPredictions,
) -> ProteinFunctionPredictions {
    ProteinFunctionPredictions {
        sift: jp.sift.and_then(decode_prediction_matrix),
        polyphen_humdiv: jp.polyphen_humdiv.and_then(decode_prediction_matrix),
        polyphen_humvar: jp.polyphen_humvar.and_then(decode_prediction_matrix),
    }
}

fn convert_exon(je: JsonExon, rank_fallback: u32) -> Exon {
    Exon {
        stable_id: je.stable_id,
        start: je.start,
        end: je.end,
        rank: je.rank.unwrap_or(rank_fallback),
        phase: je.phase,
        end_phase: je.end_phase,
    }
}

fn convert_intron(ji: JsonIntron, rank_fallback: u32) -> Intron {
    Intron {
        start: ji.start,
        end: ji.end,
        rank: ji.rank.unwrap_or(rank_fallback),
    }
}

fn attribute_value_is_true(value: &str) -> bool {
    let normalized = value.trim().to_ascii_lowercase();
    matches!(normalized.as_str(), "" | "1" | "true" | "yes")
}

/// Transcript flags are the codes of every `cds_` attribute, whatever its value:
/// VEP prints them that way (OutputFactory.pm) and its consequence predicates
/// test the attribute's presence, not its value (VariationEffect.pm
/// `_overlaps_start_codon`, `_overlaps_stop_codon`).
fn derive_flags_from_attributes(attributes: &[Attribute]) -> Vec<String> {
    let mut flags = Vec::new();
    for attr in attributes {
        if attr.code.starts_with("cds_") && !flags.contains(&attr.code) {
            flags.push(attr.code.clone());
        }
    }
    flags
}

fn derive_gencode_primary(attributes: &[Attribute]) -> bool {
    attributes.iter().any(|a| {
        (a.code == "gencode_primary" || a.code == "gencode_basic")
            && attribute_value_is_true(&a.value)
    })
}

fn convert_transcript(jt: JsonTranscript, chr_override: &str) -> Transcript {
    let chr = jt.chr.unwrap_or_else(|| chr_override.to_string());

    let exons: Vec<Exon> = jt
        .exons
        .into_iter()
        .enumerate()
        .map(|(i, e)| convert_exon(e, (i + 1) as u32))
        .collect();

    let (coding_region_start, coding_region_end) = if let Some(ref vefc) = jt.vefc {
        if let Some(ref mapper) = vefc.mapper {
            compute_coding_region(mapper)
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };

    let translation = jt.translation.map(|t| Translation {
        stable_id: t.stable_id,
        version: t.version,
        db_id: t.db_id,
        start: t.start,
        end: t.end,
        start_exon_index: 0,
        end_exon_index: exons.len().saturating_sub(1),
        seq: t.seq,
    });

    let attributes: Vec<Attribute> = jt
        .attributes
        .into_iter()
        .map(|a| Attribute {
            code: a.code,
            value: a.value,
        })
        .collect();
    let flags = derive_flags_from_attributes(&attributes);
    let gencode_primary = derive_gencode_primary(&attributes);

    let tsl = jt.tsl.or_else(|| {
        attributes.iter().find_map(|a| {
            if a.code == "TSL" {
                a.value
                    .strip_prefix("tsl")
                    .and_then(|n| n.parse::<u8>().ok())
            } else {
                None
            }
        })
    });

    // The transcript-level `introns` field prefers VEFC introns over the
    // exon-derived ones, which are only computed when needed.
    let (vefc, introns) = if let Some(v) = jt.vefc {
        let vefc_introns: Vec<Intron> = v
            .introns
            .into_iter()
            .enumerate()
            .map(|(i, intr)| convert_intron(intr, (i + 1) as u32))
            .collect();

        // `sorted_exons` is the exon list in genomic order; a cache that omits
        // it (a pruned fixture) gets the same list from `exons`.
        let sorted_exons: Vec<Exon> = if v.sorted_exons.is_empty() {
            let mut derived = exons.clone();
            derived.sort_by_key(|e| (e.start, e.end));
            for (i, e) in derived.iter_mut().enumerate() {
                e.rank = (i + 1) as u32;
            }
            derived
        } else {
            v.sorted_exons
                .into_iter()
                .enumerate()
                .map(|(i, e)| convert_exon(e, (i + 1) as u32))
                .collect()
        };

        let mapper = v.mapper.map(|m| TranscriptMapper {
            start_phase: m.start_phase,
            cdna_coding_start: m.cdna_coding_start.unwrap_or(0),
            cdna_coding_end: m.cdna_coding_end.unwrap_or(0),
            exon_coord_mapper: ExonCoordMapper::new(
                m.pairs
                    .into_iter()
                    .map(|p| MapperPair {
                        from_start: p.from_start,
                        from_end: p.from_end,
                        to_start: p.to_start,
                        to_end: p.to_end,
                        ori: p.ori,
                    })
                    .collect(),
            ),
        });

        let protein_features: Vec<ProteinFeature> = v
            .protein_features
            .into_iter()
            .map(|pf| ProteinFeature {
                start: pf.start,
                end: pf.end,
                hseqname: pf.hseqname,
                analysis: pf.analysis,
            })
            .collect();

        let tx_introns = if !vefc_introns.is_empty() {
            vefc_introns.clone()
        } else {
            derive_introns(&exons)
        };

        let protein_function_predictions = v
            .protein_function_predictions
            .map(convert_protein_function_predictions);

        let vefc = TranscriptVEFC {
            codon_table: v.codon_table.unwrap_or(1),
            five_prime_utr: v.five_prime_utr,
            three_prime_utr: v.three_prime_utr,
            translateable_seq: v.translateable_seq,
            peptide: v.peptide,
            introns: vefc_introns,
            sorted_exons,
            mapper,
            protein_features,
            protein_function_predictions,
            seq_edits: vec![],
        };

        (Some(vefc), tx_introns)
    } else {
        (None, derive_introns(&exons))
    };

    Transcript {
        stable_id: jt.stable_id,
        version: jt.version,
        db_id: jt.db_id,
        gene_stable_id: jt.gene_stable_id,
        chr,
        start: jt.start,
        end: jt.end,
        strand: jt.strand,
        biotype: jt.biotype,
        source: jt.source,
        description: jt.description,
        gene_symbol: jt.gene_symbol,
        gene_symbol_source: jt.gene_symbol_source,
        hgnc_id: jt.hgnc_id,
        gene_phenotype: jt.gene_phenotype,
        canonical: jt.canonical,
        mane_select: jt.mane_select,
        mane_plus_clinical: jt.mane_plus_clinical,
        tsl,
        appris: jt.appris,
        ccds: jt.ccds,
        protein_id: jt.protein_id,
        refseq: jt.refseq,
        swissprot: jt.swissprot,
        trembl: jt.trembl,
        uniparc: jt.uniparc,
        exons,
        introns,
        cdna_coding_start: jt.cdna_coding_start,
        cdna_coding_end: jt.cdna_coding_end,
        coding_region_start,
        coding_region_end,
        translation_start: coding_region_start,
        translation_end: coding_region_end,
        translation,
        cdna_sequence: None,
        protein_sequence: None,
        flags: flags.into(),
        gencode_primary,
        attributes,
        vefc,
        derived: Default::default(),
    }
}

/// Derive introns from a sorted list of exons.
fn derive_introns(exons: &[Exon]) -> Vec<Intron> {
    if exons.len() < 2 {
        return vec![];
    }
    let mut sorted: Vec<&Exon> = exons.iter().collect();
    sorted.sort_by_key(|e| e.start);

    sorted
        .windows(2)
        .enumerate()
        .map(|(i, pair)| Intron {
            start: pair[0].end + 1,
            end: pair[1].start - 1,
            rank: (i + 1) as u32,
        })
        .collect()
}

/// Compute coding region start/end from the mapper.
fn compute_coding_region(mapper: &JsonMapper) -> (Option<u64>, Option<u64>) {
    let cds_start = match mapper.cdna_coding_start {
        Some(v) if v > 0 => v,
        _ => return (None, None),
    };
    let cds_end = match mapper.cdna_coding_end {
        Some(v) if v > 0 => v,
        _ => return (None, None),
    };

    let mut genomic_start: Option<u64> = None;
    let mut genomic_end: Option<u64> = None;

    for pair in &mapper.pairs {
        if cds_start >= pair.from_start && cds_start <= pair.from_end {
            let offset = cds_start - pair.from_start;
            if pair.ori >= 0 {
                genomic_start = Some(pair.to_start + offset);
            } else {
                genomic_start = Some(pair.to_end - offset);
            }
        }
        if cds_end >= pair.from_start && cds_end <= pair.from_end {
            let offset = cds_end - pair.from_start;
            if pair.ori >= 0 {
                genomic_end = Some(pair.to_start + offset);
            } else {
                genomic_end = Some(pair.to_end - offset);
            }
        }
    }

    // For reverse strand, start > end in genomic coords, so swap
    if let (Some(s), Some(e)) = (genomic_start, genomic_end) {
        if s > e {
            return (Some(e), Some(s));
        }
    }

    (genomic_start, genomic_end)
}

/// Enumerate the per-chromosome JSON shard paths in a cache directory.
///
/// Directory structure: `{json_cache_dir}/transcripts/{chr}/{start}-{end}.json`.
/// This is a `read_dir` walk only: no JSON is parsed and no transcript is
/// constructed, so it is cheap enough to run eagerly at startup even when the
/// run will only ever query one chromosome.
///
/// The returned map's key set is the complete set of chromosomes the cache
/// covers. That matters beyond convenience: several Perl-parity predicates in
/// `runner.rs` / `annotator.rs` branch on whether a chromosome is *present* in
/// the cache (the `intergenic_variant` fallback, and the inter-chromosomal BND
/// mate-annotation gate), mirroring Perl VEP's own cache-loading behavior. A
/// lazy loader must therefore be able to answer "is this chromosome in the
/// cache?" without parsing it, which this function provides.
///
/// Returns an empty map when the directory is absent.
///
/// A chromosome directory holding no `.json` shards yields no key. The
/// distinction is load-bearing rather than cosmetic: `contains_key` gates the
/// `intergenic_variant` fallback and the inter-chromosomal BND mate check, so a
/// key with an empty shard list would materialize an empty index and flip both
/// predicates from "chromosome absent" (drop the row, as Perl does) to
/// "chromosome present but has no transcripts" (emit `intergenic_variant`). A
/// partially-copied cache is the realistic way to end up with an empty
/// chromosome directory.
pub fn enumerate_transcript_shards(
    json_cache_dir: &str,
) -> anyhow::Result<HashMap<String, Vec<PathBuf>>> {
    let mut by_chr: HashMap<String, Vec<PathBuf>> = HashMap::new();
    let transcripts_dir = Path::new(json_cache_dir).join("transcripts");

    if !transcripts_dir.exists() {
        info!(
            "No transcripts directory found at {}",
            transcripts_dir.display()
        );
        return Ok(by_chr);
    }

    for chr_entry in std::fs::read_dir(&transcripts_dir)
        .with_context(|| format!("Failed to read {}", transcripts_dir.display()))?
    {
        let chr_entry = chr_entry?;
        if !chr_entry.file_type()?.is_dir() {
            continue;
        }
        let chr = chr_entry.file_name().to_string_lossy().to_string();

        let mut shards: Vec<PathBuf> = Vec::new();
        for file_entry in std::fs::read_dir(chr_entry.path())? {
            let file_entry = file_entry?;
            let path = file_entry.path();
            if !is_shard_path(&path) {
                continue;
            }
            shards.push(path);
        }
        // Shards load in genomic order, whatever order the directory listing
        // returned them in, so two caches with the same content index their
        // transcripts identically.
        shards.sort_by_key(|p| shard_start(p));

        // Insert only when the directory actually holds shards, so a shard-less
        // chromosome directory reads as absent (see the note above).
        if !shards.is_empty() {
            by_chr.insert(chr, shards);
        }
    }

    Ok(by_chr)
}

/// A shard is `<start>-<end>.json`, optionally gzip-compressed as `.json.gz`.
fn is_shard_path(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".json") || name.ends_with(".json.gz")
}

/// The numeric start of a shard's range from its file name; a name that does
/// not start with a number sorts last, in name order.
fn shard_start(path: &Path) -> (u64, String) {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_string();
    let start = name
        .split('-')
        .next()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(u64::MAX);
    (start, name)
}

/// Read a shard's bytes, inflating a `.json.gz` shard.
fn read_shard(path: &Path) -> anyhow::Result<Vec<u8>> {
    let raw = std::fs::read(path).with_context(|| format!("Failed to read {}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        let mut out = Vec::new();
        std::io::Read::read_to_end(&mut flate2::read::GzDecoder::new(&raw[..]), &mut out)
            .with_context(|| format!("Failed to inflate {}", path.display()))?;
        return Ok(out);
    }
    Ok(raw)
}

/// Parse a shard's transcripts from either shard shape.
///
/// vep-cache-builder writes a bare array; a Storable-derived cache wraps the
/// transcripts in a range-keyed map (`{"<chr>": [transcript, null, ...]}`) whose
/// nulls are dropped. The first non-whitespace byte selects the shape, so a
/// shard is parsed once and a malformed shard reports the error from the parse
/// that matched its shape. Map keys are visited in sorted order so two loads of
/// the same shard yield the same transcript order.
fn parse_shard_transcripts(data: &[u8]) -> serde_json::Result<Vec<JsonTranscript>> {
    let first = data
        .iter()
        .find(|b| !matches!(b, b' ' | b'\t' | b'\n' | b'\r'))
        .copied();
    if first == Some(b'{') {
        let map: BTreeMap<String, Vec<Option<JsonTranscript>>> = serde_json::from_slice(data)?;
        return Ok(map
            .into_values()
            .flat_map(|v| v.into_iter().flatten())
            .collect());
    }
    serde_json::from_slice(data)
}

/// Parse one chromosome's JSON shards into `Transcript`s.
///
/// Shards are parsed in parallel via rayon on the global pool, not the scoped
/// pool `Runner` builds from `--fork`: both callers run outside `pool.install`
/// (the lazy path deliberately so; see [`crate::transcript_index::LazyTranscriptIndexes::prewarm`]),
/// and work only lands on a scoped pool from inside its own `install`. So
/// `--fork 1` does not make shard parsing single-threaded.
///
/// `chr` is the chromosome name the shards belong to; it is stamped onto every
/// transcript by `convert_transcript`.
pub fn load_transcripts_for_chr(chr: &str, shards: &[PathBuf]) -> anyhow::Result<Vec<Transcript>> {
    let parsed: Vec<anyhow::Result<Vec<Transcript>>> = shards
        .par_iter()
        .map(|path| -> anyhow::Result<Vec<Transcript>> {
            let data = read_shard(path)?;
            let json_transcripts = parse_shard_transcripts(&data)
                .with_context(|| format!("Failed to parse {}", path.display()))?;
            let transcripts: Vec<Transcript> = json_transcripts
                .into_iter()
                .map(|jt| convert_transcript(jt, chr))
                .collect();
            Ok(transcripts)
        })
        .collect();

    // A cache converted from VEP's region files stores a gene's transcripts in
    // every 1 Mb shard the gene overlaps, so the same transcript arrives from two
    // shards; VEP de-duplicates them by dbID on load, and so does this loader,
    // keeping the first copy. Caches without dbIDs are keyed on stable id and
    // version.
    let mut out: Vec<Transcript> = Vec::new();
    let mut seen: std::collections::HashSet<(Option<u64>, Arc<str>, Option<u32>)> =
        std::collections::HashSet::new();
    for item in parsed {
        for t in item? {
            let key = match t.db_id {
                Some(id) => (Some(id), Arc::from(""), None),
                None => (None, t.stable_id.clone(), t.version),
            };
            if seen.insert(key) {
                out.push(t);
            }
        }
    }
    Ok(out)
}

/// Load every chromosome's transcripts from the JSON cache directory.
///
/// Eager counterpart to [`enumerate_transcript_shards`] +
/// [`load_transcripts_for_chr`]. For callers that serve many requests from one
/// load (a long-lived consumer), where deferring per-chromosome parsing buys
/// nothing.
///
/// Per-chromosome counts are logged here, in sorted order, rather than inside
/// `load_transcripts_for_chr`: the loop below is parallel, so logging from the
/// closure would interleave the ~25 lines in completion order and make two runs
/// of the same cache produce different logs.
pub fn load_all_transcripts(
    json_cache_dir: &str,
) -> anyhow::Result<HashMap<String, Vec<Transcript>>> {
    let shards = enumerate_transcript_shards(json_cache_dir)?;
    // Parallel across chromosomes, not only across one chromosome's shards: a
    // serial outer loop drains the pool at every chromosome boundary and runs
    // each tail shard alone. Rayon nests safely here because this eager path
    // runs once at startup outside any `pool.install` scope; the lazy path
    // pre-warms serially for that reason (`LazyTranscriptIndexes::prewarm`).
    let by_chr: HashMap<String, Vec<Transcript>> = shards
        .par_iter()
        .map(
            |(chr, paths)| -> anyhow::Result<(String, Vec<Transcript>)> {
                Ok((chr.clone(), load_transcripts_for_chr(chr, paths)?))
            },
        )
        .collect::<anyhow::Result<HashMap<_, _>>>()?;
    info!("Loaded transcripts for {} chromosomes", by_chr.len());
    let mut chrs: Vec<&String> = by_chr.keys().collect();
    chrs.sort();
    for chr in chrs {
        info!("  chr{}: {} transcripts", chr, by_chr[chr].len());
    }
    Ok(by_chr)
}

/// Deserialize helper for fields that can be "", "0", 0, or 1 in the variation JSON.
fn deserialize_string_or_num_as_u8<'de, D>(deserializer: D) -> Result<u8, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum StringOrNum {
        Num(i64),
        Float(f64),
        Str(String),
    }
    let v: Option<StringOrNum> = Option::deserialize(deserializer)?;
    match v {
        None => Ok(0),
        Some(StringOrNum::Num(n)) => Ok(n as u8),
        Some(StringOrNum::Float(f)) => Ok(f as u8),
        Some(StringOrNum::Str(s)) if s.is_empty() => Ok(0),
        Some(StringOrNum::Str(s)) => s.parse::<u8>().or(Ok(0)),
    }
}

fn deserialize_var_start<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Val {
        Num(u64),
        Float(f64),
        Str(String),
    }
    match Val::deserialize(deserializer)? {
        Val::Num(n) => Ok(n),
        Val::Float(f) => Ok(f as u64),
        Val::Str(s) => s.parse::<u64>().map_err(serde::de::Error::custom),
    }
}

fn deserialize_optional_var_u64<'de, D>(deserializer: D) -> Result<Option<u64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Val {
        Num(u64),
        Float(f64),
        Str(String),
    }
    let v: Option<Val> = Option::deserialize(deserializer)?;
    match v {
        None => Ok(None),
        Some(Val::Str(s)) if s.is_empty() => Ok(None),
        Some(Val::Str(s)) => s.parse::<u64>().map(Some).map_err(serde::de::Error::custom),
        Some(Val::Num(n)) => Ok(Some(n)),
        Some(Val::Float(f)) => Ok(Some(f as u64)),
    }
}

fn deserialize_optional_var_f64<'de, D>(deserializer: D) -> Result<Option<f64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Val {
        Num(f64),
        Str(String),
    }
    let v: Option<Val> = Option::deserialize(deserializer)?;
    match v {
        None => Ok(None),
        Some(Val::Str(s)) if s.is_empty() => Ok(None),
        Some(Val::Str(s)) => s.parse::<f64>().map(Some).map_err(serde::de::Error::custom),
        Some(Val::Num(n)) => Ok(Some(n)),
    }
}

fn deserialize_var_strand<'de, D>(deserializer: D) -> Result<i8, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Val {
        Num(i64),
        Float(f64),
        Str(String),
    }
    let v: Option<Val> = Option::deserialize(deserializer)?;
    match v {
        None => Ok(1),
        Some(Val::Num(n)) => Ok(n as i8),
        Some(Val::Float(f)) => Ok(f as i8),
        Some(Val::Str(s)) if s.is_empty() => Ok(1),
        Some(Val::Str(s)) => s.parse::<i8>().or(Ok(1)),
    }
}

/// JSON representation of a variation record from the converted cache.
#[derive(Deserialize)]
#[allow(non_snake_case)] // The cache carries Perl's mixed-case field names (e.g. "gnomADe_AFR")
struct JsonVariation {
    #[serde(default)]
    variation_name: Option<String>,
    #[serde(default, deserialize_with = "deserialize_string_or_num_as_u8")]
    failed: u8,
    #[serde(default, deserialize_with = "deserialize_string_or_num_as_u8")]
    somatic: u8,
    #[serde(deserialize_with = "deserialize_var_start")]
    start: u64,
    #[serde(default, deserialize_with = "deserialize_optional_var_u64")]
    end: Option<u64>,
    #[serde(default)]
    allele_string: Option<String>,
    #[serde(default, deserialize_with = "deserialize_var_strand")]
    strand: i8,
    #[serde(default)]
    minor_allele: Option<String>,
    #[serde(default, deserialize_with = "deserialize_optional_var_f64")]
    minor_allele_freq: Option<f64>,
    #[serde(default)]
    clin_sig: Option<String>,
    #[serde(default, deserialize_with = "deserialize_string_or_num_as_u8")]
    phenotype_or_disease: u8,
    #[serde(default)]
    pubmed: Option<String>,
    #[serde(default)]
    var_synonyms: Option<String>,
    #[serde(default, rename = "AF")]
    af: Option<String>,
    #[serde(default, rename = "AFR")]
    afr: Option<String>,
    #[serde(default, rename = "AMR")]
    amr: Option<String>,
    #[serde(default, rename = "EAS")]
    eas: Option<String>,
    #[serde(default, rename = "EUR")]
    eur: Option<String>,
    #[serde(default, rename = "SAS")]
    sas: Option<String>,
    #[serde(default, rename = "AA")]
    aa: Option<String>,
    #[serde(default, rename = "EA")]
    ea: Option<String>,
    #[serde(default)]
    gnomADe: Option<String>,
    #[serde(default)]
    gnomADe_AFR: Option<String>,
    #[serde(default)]
    gnomADe_AMR: Option<String>,
    #[serde(default)]
    gnomADe_ASJ: Option<String>,
    #[serde(default)]
    gnomADe_EAS: Option<String>,
    #[serde(default)]
    gnomADe_FIN: Option<String>,
    #[serde(default)]
    gnomADe_MID: Option<String>,
    #[serde(default)]
    gnomADe_NFE: Option<String>,
    #[serde(default)]
    gnomADe_REMAINING: Option<String>,
    #[serde(default)]
    gnomADe_SAS: Option<String>,
    #[serde(default)]
    gnomADg: Option<String>,
    #[serde(default)]
    gnomADg_AFR: Option<String>,
    #[serde(default)]
    gnomADg_AMR: Option<String>,
    #[serde(default)]
    gnomADg_ASJ: Option<String>,
    #[serde(default)]
    gnomADg_EAS: Option<String>,
    #[serde(default)]
    gnomADg_FIN: Option<String>,
    #[serde(default)]
    gnomADg_MID: Option<String>,
    #[serde(default)]
    gnomADg_NFE: Option<String>,
    #[serde(default)]
    gnomADg_REMAINING: Option<String>,
    #[serde(default)]
    gnomADg_SAS: Option<String>,
}

fn convert_variation(jv: JsonVariation) -> CachedVariation {
    let mut frequencies = HashMap::new();

    macro_rules! insert_freq {
        ($key:expr, $val:expr) => {
            if let Some(ref v) = $val {
                if !v.is_empty() {
                    frequencies.insert($key.to_string(), v.clone());
                }
            }
        };
    }

    insert_freq!("AF", jv.af);
    insert_freq!("AFR", jv.afr);
    insert_freq!("AMR", jv.amr);
    insert_freq!("EAS", jv.eas);
    insert_freq!("EUR", jv.eur);
    insert_freq!("SAS", jv.sas);
    insert_freq!("AA", jv.aa);
    insert_freq!("EA", jv.ea);
    insert_freq!("gnomADe", jv.gnomADe);
    insert_freq!("gnomADe_AFR", jv.gnomADe_AFR);
    insert_freq!("gnomADe_AMR", jv.gnomADe_AMR);
    insert_freq!("gnomADe_ASJ", jv.gnomADe_ASJ);
    insert_freq!("gnomADe_EAS", jv.gnomADe_EAS);
    insert_freq!("gnomADe_FIN", jv.gnomADe_FIN);
    insert_freq!("gnomADe_MID", jv.gnomADe_MID);
    insert_freq!("gnomADe_NFE", jv.gnomADe_NFE);
    insert_freq!("gnomADe_REMAINING", jv.gnomADe_REMAINING);
    insert_freq!("gnomADe_SAS", jv.gnomADe_SAS);
    insert_freq!("gnomADg", jv.gnomADg);
    insert_freq!("gnomADg_AFR", jv.gnomADg_AFR);
    insert_freq!("gnomADg_AMR", jv.gnomADg_AMR);
    insert_freq!("gnomADg_ASJ", jv.gnomADg_ASJ);
    insert_freq!("gnomADg_EAS", jv.gnomADg_EAS);
    insert_freq!("gnomADg_FIN", jv.gnomADg_FIN);
    insert_freq!("gnomADg_MID", jv.gnomADg_MID);
    insert_freq!("gnomADg_NFE", jv.gnomADg_NFE);
    insert_freq!("gnomADg_REMAINING", jv.gnomADg_REMAINING);
    insert_freq!("gnomADg_SAS", jv.gnomADg_SAS);

    let start = jv.start;
    let end = jv.end.unwrap_or(start);

    let clin_sig = jv.clin_sig.filter(|s| !s.is_empty());
    let pubmed = jv.pubmed.filter(|s| !s.is_empty());
    let minor_allele = jv.minor_allele.filter(|s| !s.is_empty());
    let var_synonyms = jv.var_synonyms.filter(|s| !s.is_empty());

    CachedVariation {
        variation_name: jv.variation_name.unwrap_or_default(),
        failed: jv.failed,
        somatic: jv.somatic,
        start,
        end,
        allele_string: jv.allele_string.unwrap_or_default(),
        strand: jv.strand,
        minor_allele,
        minor_allele_freq: jv.minor_allele_freq,
        clin_sig,
        phenotype_or_disease: jv.phenotype_or_disease,
        pubmed,
        frequencies,
        var_synonyms,
    }
}

/// Container for loaded variation data with a position index for fast lookup.
///
/// The `position_index` maps `(chromosome, position) -> Vec<index>` into the
/// per-chromosome `variations` vec, enabling O(1) position-based lookups in
/// the variation matching pipeline.
pub struct VariationData {
    /// Variations keyed by chromosome.
    pub variations: HashMap<String, Vec<CachedVariation>>,
    /// Position index: chromosome -> (position -> vec of indices into variations vec).
    pub position_index: HashMap<String, FxHashMap<u64, Vec<usize>>>,
}

/// Load all variations from the JSON cache directory.
///
/// Directory structure: `{json_cache_dir}/variations/{chr}/{start}-{end}.json`
///
/// JSON files are parsed in parallel using the rayon thread pool (same
/// strategy as `load_all_transcripts`).
pub fn load_all_variations(json_cache_dir: &str) -> anyhow::Result<VariationData> {
    let mut variations: HashMap<String, Vec<CachedVariation>> = HashMap::new();
    let mut position_index: HashMap<String, FxHashMap<u64, Vec<usize>>> = HashMap::new();
    let variations_dir = Path::new(json_cache_dir).join("variations");

    if !variations_dir.exists() {
        info!(
            "No variations directory found at {}",
            variations_dir.display()
        );
        return Ok(VariationData {
            variations,
            position_index,
        });
    }

    let mut files: Vec<(String, PathBuf)> = Vec::new();
    for chr_entry in std::fs::read_dir(&variations_dir)
        .with_context(|| format!("Failed to read {}", variations_dir.display()))?
    {
        let chr_entry = chr_entry?;
        if !chr_entry.file_type()?.is_dir() {
            continue;
        }
        let chr = chr_entry.file_name().to_string_lossy().to_string();

        for file_entry in std::fs::read_dir(chr_entry.path())? {
            let file_entry = file_entry?;
            let path = file_entry.path();
            if path.extension().is_none_or(|e| e != "json") {
                continue;
            }
            files.push((chr.clone(), path));
        }
    }

    let parsed: Vec<anyhow::Result<(String, Vec<CachedVariation>)>> = files
        .par_iter()
        .map(
            |(chr, path)| -> anyhow::Result<(String, Vec<CachedVariation>)> {
                let data = std::fs::read(path)
                    .with_context(|| format!("Failed to read {}", path.display()))?;
                let json_variations: Vec<JsonVariation> = serde_json::from_slice(&data)
                    .with_context(|| format!("Failed to parse {}", path.display()))?;
                let converted: Vec<CachedVariation> =
                    json_variations.into_iter().map(convert_variation).collect();
                Ok((chr.clone(), converted))
            },
        )
        .collect();

    for item in parsed {
        let (chr, converted) = item?;
        let chr_vars = variations.entry(chr.clone()).or_default();
        let chr_idx = position_index.entry(chr).or_default();
        let base_idx = chr_vars.len();

        for (i, cv) in converted.into_iter().enumerate() {
            let pos = cv.start;
            chr_idx.entry(pos).or_default().push(base_idx + i);
            chr_vars.push(cv);
        }
    }

    info!("Loaded variations for {} chromosomes", variations.len());
    for (chr, vars) in &variations {
        info!("  chr{}: {} variations", chr, vars.len());
    }

    Ok(VariationData {
        variations,
        position_index,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_cache_missing_dir() {
        let result = load_all_transcripts("/nonexistent/path");
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    /// The parallel outer loop must produce the same map a serial one would.
    ///
    /// `load_all_transcripts` fans out across chromosomes with `par_iter` and
    /// `load_transcripts_for_chr` fans out across each chromosome's shards, so a
    /// chromosome's transcripts arrive from several workers. This asserts nothing is
    /// dropped, duplicated, or attributed to the wrong chromosome; multiple shards
    /// per chromosome and differing counts make a mis-keyed result detectable.
    #[test]
    fn test_load_all_transcripts_parallel_matches_expected_grouping() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_parallel_")
            .tempdir()
            .unwrap();
        // chr21: 3 shards x 2 transcripts; chr1: 2 shards x 1; chr2: 1 shard x 4.
        for (chr, shards, per_shard) in [("21", 3usize, 2usize), ("1", 2, 1), ("2", 1, 4)] {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            for shard in 0..shards {
                let entries: Vec<String> = (0..per_shard)
                    .map(|i| {
                        format!(
                            r#"{{"stable_id":"ENST{chr}_{shard}_{i}","gene_stable_id":"ENSG{chr}","chr":"{chr}","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#,
                            start = 100 + (shard * 1000 + i * 10) as u64,
                            end = 200 + (shard * 1000 + i * 10) as u64
                        )
                    })
                    .collect();
                std::fs::write(
                    chr_dir.join(format!("{shard}-shard.json")),
                    format!("[{}]", entries.join(",")),
                )
                .unwrap();
            }
        }

        let by_chr = load_all_transcripts(dir.path().to_str().unwrap()).unwrap();

        let mut counts: Vec<(String, usize)> = by_chr
            .iter()
            .map(|(chr, txs)| (chr.clone(), txs.len()))
            .collect();
        counts.sort();
        assert_eq!(
            counts,
            vec![
                ("1".to_string(), 2),
                ("2".to_string(), 4),
                ("21".to_string(), 6),
            ],
            "every shard's transcripts must land under its own chromosome exactly once"
        );

        // Every transcript must carry the chromosome it was filed under: the `chr`
        // stamp is applied inside the parallel closure, so a mixed-up pairing would
        // preserve the counts above while corrupting the data.
        for (chr, txs) in &by_chr {
            assert!(
                txs.iter().all(|t| &*t.chr == chr.as_str()),
                "chr{chr} transcripts must all be stamped chr{chr}"
            );
            let mut ids: Vec<&str> = txs.iter().map(|t| &*t.stable_id).collect();
            ids.sort();
            let unique = ids.len();
            ids.dedup();
            assert_eq!(unique, ids.len(), "chr{chr} must have no duplicate ids");
        }
    }

    /// A cache converted from VEP's region files carries a gene's transcripts in
    /// every shard the gene overlaps. The loader must keep one copy, and the copy
    /// it keeps must not depend on directory listing order.
    #[test]
    fn a_transcript_present_in_two_shards_loads_once() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_dedup_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let straddler = r#"{"stable_id":"ENST_STRADDLE","dbID":"77","version":3,"gene_stable_id":"ENSG_S","chr":"21","start":990000,"end":1010000,"strand":1,"biotype":"protein_coding","source":"ensembl"}"#;
        let only_a = r#"{"stable_id":"ENST_A","dbID":"1","gene_stable_id":"ENSG_A","chr":"21","start":100,"end":900,"strand":1,"biotype":"protein_coding","source":"ensembl"}"#;
        let no_dbid_1 = r#"{"stable_id":"ENST_NODB","version":1,"gene_stable_id":"ENSG_N","chr":"21","start":1500000,"end":1500900,"strand":1,"biotype":"protein_coding","source":"ensembl"}"#;
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!("[{only_a},{straddler}]"),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("1000001-2000000.json"),
            format!("[{straddler},{no_dbid_1}]"),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("2000001-3000000.json"),
            format!("[{no_dbid_1}]"),
        )
        .unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let names: Vec<String> = shards["21"]
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        assert_eq!(
            names,
            vec![
                "1-1000000.json",
                "1000001-2000000.json",
                "2000001-3000000.json"
            ]
        );

        let txs = load_transcripts_for_chr("21", &shards["21"]).unwrap();
        let mut ids: Vec<&str> = txs.iter().map(|t| &*t.stable_id).collect();
        ids.sort();
        assert_eq!(ids, vec!["ENST_A", "ENST_NODB", "ENST_STRADDLE"]);
    }

    /// A shard may be stored gzip-compressed as `.json.gz`; it loads like a plain shard.
    #[test]
    fn gzipped_shards_load_like_plain_ones() {
        use std::io::Write as _;
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_gz_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let body = r#"[{"stable_id":"ENST_GZ","dbID":"5","gene_stable_id":"ENSG_G","chr":"21","start":100,"end":900,"strand":1,"biotype":"protein_coding","source":"ensembl"}]"#;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        std::fs::write(chr_dir.join("1-1000000.json.gz"), enc.finish().unwrap()).unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        assert_eq!(shards["21"].len(), 1);
        let txs = load_transcripts_for_chr("21", &shards["21"]).unwrap();
        assert_eq!(txs.len(), 1);
        assert_eq!(&*txs[0].stable_id, "ENST_GZ");
    }

    /// A Storable-derived shard wraps its transcripts in a range-keyed map with
    /// null placeholders; it loads beside an array shard, nulls dropped, keys in
    /// sorted order.
    #[test]
    fn map_shaped_shards_load_beside_array_shards() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_map_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let tx = |id: &str, dbid: u32| {
            format!(
                r#"{{"stable_id":"{id}","dbID":{dbid},"gene_stable_id":"ENSG_M","start":100,"end":900,"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#
            )
        };
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!("[{}]", tx("ENST_ARR", 1)),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("1000001-2000000.json"),
            format!(
                " \n{{\"b\": [null, {}, null], \"a\": [{}, null]}}",
                tx("ENST_MAP_B", 2),
                tx("ENST_MAP_A", 3)
            ),
        )
        .unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let txs = load_transcripts_for_chr("21", &shards["21"]).unwrap();
        let ids: Vec<&str> = txs.iter().map(|t| &*t.stable_id).collect();
        assert_eq!(ids, vec!["ENST_ARR", "ENST_MAP_A", "ENST_MAP_B"]);
        assert!(txs.iter().all(|t| t.chr == "21"));
    }

    /// A truncated map shard fails naming the shard and the map parse's own
    /// error, not an "expected a sequence" complaint from an array attempt.
    #[test]
    fn malformed_map_shard_reports_the_map_parse_error() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_badmap_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            r#"{"21": [{"stable_id": "ENST_TRUNC", "gene_stable_id": "#,
        )
        .unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let err = load_transcripts_for_chr("21", &shards["21"]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("1-1000000.json"), "{msg}");
        assert!(msg.contains("EOF"), "{msg}");
        assert!(!msg.contains("expected a sequence"), "{msg}");
    }

    #[test]
    fn test_load_variations_missing_dir() {
        let result = load_all_variations("/nonexistent/path");
        assert!(result.is_ok());
        let data = result.unwrap();
        assert!(data.variations.is_empty());
        assert!(data.position_index.is_empty());
    }

    #[test]
    fn test_convert_variation_basic() {
        let json = r#"{
            "variation_name": "rs123",
            "start": "100",
            "allele_string": "A/G",
            "strand": 1,
            "failed": "",
            "somatic": "",
            "phenotype_or_disease": 0,
            "end": "",
            "AFR": "G:0.003",
            "gnomADe": "G:0.0005"
        }"#;
        let jv: JsonVariation = serde_json::from_str(json).unwrap();
        let cv = convert_variation(jv);
        assert_eq!(cv.variation_name, "rs123");
        assert_eq!(cv.start, 100);
        assert_eq!(cv.allele_string, "A/G");
        assert_eq!(cv.strand, 1);
        assert_eq!(cv.failed, 0);
        assert_eq!(cv.somatic, 0);
        assert_eq!(cv.frequencies.get("AFR").unwrap(), "G:0.003");
        assert_eq!(cv.frequencies.get("gnomADe").unwrap(), "G:0.0005");
    }

    #[test]
    fn test_convert_variation_empty_optionals() {
        let json = r#"{
            "variation_name": "rs456",
            "start": "200",
            "allele_string": "C/T",
            "strand": "",
            "failed": "",
            "somatic": "",
            "phenotype_or_disease": "",
            "end": "",
            "clin_sig": "",
            "pubmed": "",
            "minor_allele": "",
            "minor_allele_freq": ""
        }"#;
        let jv: JsonVariation = serde_json::from_str(json).unwrap();
        let cv = convert_variation(jv);
        assert_eq!(cv.strand, 1); // default for empty
        assert!(cv.clin_sig.is_none());
        assert!(cv.pubmed.is_none());
        assert!(cv.minor_allele.is_none());
        assert!(cv.minor_allele_freq.is_none());
    }

    /// The numeric helpers accept a JSON number, a float, a string holding a
    /// number and (where optional) an empty string or null, and reject the rest.
    #[test]
    fn numeric_helpers_accept_number_string_float_and_reject_junk() {
        macro_rules! de {
            ($f:ident, $json:expr) => {{
                let mut d = serde_json::Deserializer::from_str($json);
                $f(&mut d)
            }};
        }

        assert_eq!(de!(deserialize_string_or_number, "17").unwrap(), 17);
        assert_eq!(de!(deserialize_string_or_number, "\"17\"").unwrap(), 17);
        assert_eq!(de!(deserialize_string_or_number, "17.9").unwrap(), 17);
        assert_eq!(de!(deserialize_string_or_number, "-4").unwrap(), 0);
        assert!(de!(deserialize_string_or_number, "\"\"").is_err());
        assert!(de!(deserialize_string_or_number, "\"x\"").is_err());
        assert!(de!(deserialize_string_or_number, "true").is_err());
        assert!(de!(deserialize_string_or_number, "null").is_err());

        assert_eq!(
            de!(deserialize_optional_string_or_number, "279463").unwrap(),
            Some(279463)
        );
        assert_eq!(
            de!(deserialize_optional_string_or_number, "\"279463\"").unwrap(),
            Some(279463)
        );
        assert_eq!(
            de!(deserialize_optional_string_or_number, "\"\"").unwrap(),
            None
        );
        assert_eq!(
            de!(deserialize_optional_string_or_number, "null").unwrap(),
            None
        );
        assert_eq!(
            de!(deserialize_optional_string_or_number, "-1").unwrap(),
            Some(0)
        );
        assert!(de!(deserialize_optional_string_or_number, "\"x\"").is_err());

        assert_eq!(de!(deserialize_optional_u32, "13").unwrap(), Some(13));
        assert_eq!(de!(deserialize_optional_u32, "\"13\"").unwrap(), Some(13));
        assert_eq!(de!(deserialize_optional_u32, "\"\"").unwrap(), None);
        assert_eq!(de!(deserialize_optional_u32, "null").unwrap(), None);
        assert!(de!(deserialize_optional_u32, "true").is_err());

        assert_eq!(de!(deserialize_optional_u8, "2").unwrap(), Some(2));
        assert_eq!(de!(deserialize_optional_u8, "\"2\"").unwrap(), Some(2));
        assert_eq!(de!(deserialize_optional_u8, "true").unwrap(), Some(1));
        assert_eq!(de!(deserialize_optional_u8, "false").unwrap(), Some(0));
        assert_eq!(de!(deserialize_optional_u8, "\"\"").unwrap(), None);
        assert_eq!(de!(deserialize_optional_u8, "null").unwrap(), None);

        assert_eq!(de!(deserialize_i8_from_number, "-1").unwrap(), -1);
        assert_eq!(de!(deserialize_i8_from_number, "2").unwrap(), 2);
        assert_eq!(de!(deserialize_i8_from_number, "\"-1\"").unwrap(), -1);
        assert_eq!(de!(deserialize_i8_from_number, "\"x\"").unwrap(), 0);
        assert!(de!(deserialize_i8_from_number, "1.5").is_err());
        assert!(de!(deserialize_i8_from_number, "null").is_err());

        assert_eq!(de!(deserialize_strand, "-1").unwrap(), Strand::Reverse);
        assert_eq!(de!(deserialize_strand, "1").unwrap(), Strand::Forward);
        assert_eq!(de!(deserialize_strand, "\"-1\"").unwrap(), Strand::Reverse);
        assert_eq!(de!(deserialize_strand, "\"x\"").unwrap(), Strand::Forward);

        assert!(de!(deserialize_bool_from_int, "1").unwrap());
        assert!(!de!(deserialize_bool_from_int, "0").unwrap());
        assert!(de!(deserialize_bool_from_int, "true").unwrap());
        assert!(de!(deserialize_bool_from_int, "\"1\"").is_err());
    }

    /// The display path reads the transcript intron list and the consequence
    /// path the VEFC list; the loader fills the first from the second, so the
    /// two frameshift-intron facts agree on every loaded transcript.
    #[test]
    fn loaded_transcript_frameshift_facts_agree() {
        let with_vefc = r#"{"stable_id":"ENST1","gene_stable_id":"ENSG1","start":100,"end":300,
            "strand":1,"biotype":"protein_coding","source":"ensembl",
            "exons":[{"start":100,"end":200},{"start":213,"end":300}],
            "variation_effect_feature_cache":{"introns":[{"start":201,"end":212}]}}"#;
        let jt: JsonTranscript = serde_json::from_str(with_vefc).unwrap();
        let facts = *convert_transcript(jt, "1").facts();
        assert!(facts.has_frameshift_intron);
        assert_eq!(
            facts.has_frameshift_intron,
            facts.vefc_has_frameshift_intron
        );

        let without_vefc = r#"{"stable_id":"ENST2","gene_stable_id":"ENSG2","start":100,"end":300,
            "strand":-1,"biotype":"protein_coding","source":"ensembl",
            "exons":[{"start":215,"end":300},{"start":100,"end":200}]}"#;
        let jt: JsonTranscript = serde_json::from_str(without_vefc).unwrap();
        let facts = *convert_transcript(jt, "1").facts();
        assert!(!facts.has_frameshift_intron);
        assert_eq!(
            facts.has_frameshift_intron,
            facts.vefc_has_frameshift_intron
        );
    }

    #[test]
    fn test_derive_flags_from_attributes() {
        let attrs = vec![
            Attribute {
                code: "cds_start_NF".into(),
                value: "1".into(),
            },
            Attribute {
                code: "cds_end_NF".into(),
                value: "0".into(),
            },
            Attribute {
                code: "mRNA_end_NF".into(),
                value: "true".into(),
            },
        ];
        let flags = derive_flags_from_attributes(&attrs);
        assert_eq!(
            flags,
            vec!["cds_start_NF".to_string(), "cds_end_NF".to_string()]
        );
    }

    #[test]
    fn test_derive_gencode_primary_from_attributes() {
        let attrs = vec![Attribute {
            code: "gencode_basic".into(),
            value: "1".into(),
        }];
        assert!(derive_gencode_primary(&attrs));

        let attrs_false = vec![Attribute {
            code: "gencode_primary".into(),
            value: "0".into(),
        }];
        assert!(!derive_gencode_primary(&attrs_false));
    }
}
