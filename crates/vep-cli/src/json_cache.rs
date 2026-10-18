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

use vep_core::cache_info::CacheInfo;
use vep_core::coordinate::Strand;
use vep_core::transcript::{
    Attribute, Exon, ExonCoordMapper, Intron, MapperPair, PredictionMatrix, ProteinFeature,
    ProteinFunctionPredictions, SeqEdit, Transcript, TranscriptMapper, TranscriptVEFC, Translation,
};
use vep_core::variation::{CachedVariation, ClinVarDetail};

use crate::config::CliError;

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
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    gene_version: Option<u32>,
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
    uniprot_isoform: Option<String>,
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
    #[serde(default)]
    seq_edits: Vec<JsonSeqEdit>,
}

/// A translation SeqEdit as the Storable cache stores it (`Bio::EnsEMBL::SeqEdit`
/// over `Translation::get_all_SeqEdits`): `start` and `end` are protein
/// positions, `alt_seq` the residues written over them. A cache without the key
/// loads none, and the reference peptide is then read from the cached
/// translation.
#[derive(Deserialize)]
struct JsonSeqEdit {
    #[serde(deserialize_with = "deserialize_string_or_number")]
    start: u64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(default)]
    alt_seq: String,
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

/// Deserialize `analysis`, the protein feature's display label: a string, or
/// the Storable `Bio::EnsEMBL::Analysis` object the dumper reduces to its one
/// `_display_label` key (`display_label` is the same key as the API's accessor
/// spells it).
fn deserialize_analysis<'de, D>(deserializer: D) -> Result<Option<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de;

    struct AnalysisVisitor;
    impl<'de> de::Visitor<'de> for AnalysisVisitor {
        type Value = Option<String>;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a string, null, or an object with display_label or _display_label")
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
                if key == "display_label" || key == "_display_label" {
                    display_label = val.as_str().filter(|s| !s.is_empty()).map(str::to_string);
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

/// The value of the first attribute carrying `code`: VEP reads MANE, TSL and
/// APPRIS from the attribute list with a `grep` on the code that keeps its
/// first match (OutputFactory.pm `BaseTranscriptVariationAllele_to_output_hash`),
/// so a later duplicate is never read.
fn attribute_value<'a>(attributes: &'a [Attribute], code: &str) -> Option<&'a str> {
    attributes
        .iter()
        .find(|a| a.code == code)
        .map(|a| a.value.as_str())
}

/// Whether any attribute carries `code`, whatever its value.
fn has_attribute(attributes: &[Attribute], code: &str) -> bool {
    attributes.iter().any(|a| a.code == code)
}

/// The support level a `TSL` attribute's value carries: VEP's `m/tsl(\d+)/`
/// followed by `if $1`, so `tsl1 (assigned to previous version 8)` is 1, and
/// `tslNA` and a zero are none. The first `tsl` followed by a digit is the
/// match, as the regex's is; a zero there is none and no later run is read.
fn tsl_from_attribute_value(value: &str) -> Option<u8> {
    let digits = value.match_indices("tsl").find_map(|(at, tsl)| {
        let rest = &value[at + tsl.len()..];
        let digit_run_len = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        (digit_run_len > 0).then_some(&rest[..digit_run_len])
    })?;
    digits.parse::<u8>().ok().filter(|&level| level != 0)
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

    // GENCODE set membership is the presence of the attribute, whatever its value:
    // VEP's transcript filter and its `GENCODE_PRIMARY` column test
    // `grep {$_->{code} eq 'gencode_basic'}` and `'gencode_primary'` alone
    // (AnnotationType/Transcript.pm `filter_transcript`), and the dumped values
    // are the labels `GENCODE basic` and `GENCODE primary`.
    let gencode_basic = has_attribute(&attributes, "gencode_basic");
    let gencode_primary = has_attribute(&attributes, "gencode_primary");

    // A Storable-derived cache carries MANE, TSL and APPRIS only as transcript
    // attributes; a top-level key, when a cache writes one, is read first.
    let tsl = jt
        .tsl
        .or_else(|| attribute_value(&attributes, "TSL").and_then(tsl_from_attribute_value));
    let mane_select = jt
        .mane_select
        .or_else(|| attribute_value(&attributes, "MANE_Select").map(str::to_string));
    let mane_plus_clinical = jt
        .mane_plus_clinical
        .or_else(|| attribute_value(&attributes, "MANE_Plus_Clinical").map(str::to_string));
    let appris = jt
        .appris
        .or_else(|| attribute_value(&attributes, "appris").map(str::to_string));

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

        let seq_edits: Vec<SeqEdit> = v
            .seq_edits
            .into_iter()
            .map(|e| SeqEdit {
                start: e.start,
                end: e.end,
                alt_seq: e.alt_seq,
            })
            .collect();

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
            seq_edits,
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
        gene_version: jt.gene_version,
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
        mane_select,
        mane_plus_clinical,
        tsl,
        appris,
        ccds: jt.ccds,
        protein_id: jt.protein_id,
        refseq: jt.refseq,
        swissprot: jt.swissprot,
        trembl: jt.trembl,
        uniparc: jt.uniparc,
        uniprot_isoform: jt.uniprot_isoform,
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
        gencode_basic,
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
/// Returns an empty map when the directory is absent. A `transcripts`
/// directory under which no chromosome holds a shard is refused under the error
/// contract, naming the cache's serialiser: a Sereal cache converts to no
/// shard at all, and a copy that landed no shard reads the same.
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

    if by_chr.is_empty() {
        return Err(CliError::new(
            "cache",
            format!(
                "no transcript shards under {} (serialiser_type {})",
                transcripts_dir.display(),
                cache_serialiser_type(json_cache_dir)
            ),
        )
        .into());
    }

    Ok(by_chr)
}

/// The `serialiser_type` the cache's `info.json` declares; a cache that
/// declares none was written by Storable.
fn cache_serialiser_type(json_cache_dir: &str) -> String {
    std::fs::read_to_string(Path::new(json_cache_dir).join("info.json"))
        .ok()
        .and_then(|text| serde_json::from_str::<CacheInfo>(&text).ok())
        .and_then(|info| info.serialiser_type)
        .unwrap_or_else(|| "storable".to_string())
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
    clin_sig_allele: Option<String>,
    #[serde(default)]
    clin_sig_ref_allele: Option<String>,
    #[serde(default)]
    clinical_impact: Option<String>,
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
    gnomADg_AMI: Option<String>,
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
    insert_freq!("gnomADg_AMI", jv.gnomADg_AMI);
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
    let clinvar = ClinVarDetail {
        clin_sig_allele: jv.clin_sig_allele.filter(|s| !s.is_empty()),
        clin_sig_ref_allele: jv.clin_sig_ref_allele.filter(|s| !s.is_empty()),
        clinical_impact: jv.clinical_impact.filter(|s| !s.is_empty()),
    };
    let clinvar = (clinvar != ClinVarDetail::default()).then(|| Box::new(clinvar));
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
        clinvar,
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
    /// every shard the gene overlaps. The loader keeps one copy per `dbID`, as
    /// VEP's cache reader does, and keys a transcript without one on its stable
    /// id and version; the copy it keeps does not depend on directory listing
    /// order.
    #[test]
    fn dedup_by_dbid_then_stable_id_version() {
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

    /// A transcript's `seq_edits` load in protein coordinates with the residues
    /// written over them; a cache written without the key loads with none and keeps
    /// the cached-translation reading of the reference peptide.
    #[test]
    fn seq_edits_loaded() {
        let with_edits = r#"{"stable_id":"ENST1","gene_stable_id":"ENSG1","start":100,"end":300,
            "strand":1,"biotype":"protein_coding","source":"ensembl",
            "exons":[{"start":100,"end":300}],
            "variation_effect_feature_cache":{"peptide":"MAU","translateable_seq":"CTGGCTTGA",
            "seq_edits":[{"code":"initial_met","name":"Initial methionine","start":1,"end":1,"alt_seq":"M"},
                         {"code":"_selenocysteine","start":"3","end":"3","alt_seq":"U"}]}}"#;
        let jt: JsonTranscript = serde_json::from_str(with_edits).unwrap();
        let vefc = convert_transcript(jt, "1").vefc.unwrap();
        let edits: Vec<(u64, u64, &str)> = vefc
            .seq_edits
            .iter()
            .map(|e| (e.start, e.end, e.alt_seq.as_str()))
            .collect();
        assert_eq!(edits, vec![(1, 1, "M"), (3, 3, "U")]);

        let without_key = r#"{"stable_id":"ENST2","gene_stable_id":"ENSG2","start":100,"end":300,
            "strand":1,"biotype":"protein_coding","source":"ensembl",
            "exons":[{"start":100,"end":300}],
            "variation_effect_feature_cache":{"peptide":"MAU","translateable_seq":"CTGGCTTGA"}}"#;
        let jt: JsonTranscript = serde_json::from_str(without_key).unwrap();
        assert!(convert_transcript(jt, "1")
            .vefc
            .unwrap()
            .seq_edits
            .is_empty());
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

    /// A transcript parsed from one JSON object, filed under chromosome `21`.
    fn transcript_from(json: &str) -> Transcript {
        let jt: JsonTranscript = serde_json::from_str(json).unwrap_or_else(|e| panic!("{e}"));
        convert_transcript(jt, "21")
    }

    /// The smallest transcript object the loader accepts, with `extra` JSON
    /// members spliced in after `source`.
    fn minimal(extra: &str) -> String {
        format!(
            r#"{{"stable_id":"ENST00000400375","gene_stable_id":"ENSG00000215455","start":44538981,"end":44540195,"strand":-1,"biotype":"protein_coding","source":"havana"{extra}}}"#
        )
    }

    /// One `attributes` member from `(code, value)` pairs, each as the cache
    /// writes it (with the `name` the loader ignores).
    fn attributes(pairs: &[(&str, &str)]) -> String {
        let items: Vec<String> = pairs
            .iter()
            .map(|(code, value)| {
                format!(r#"{{"code":"{code}","name":"{code} name","value":"{value}"}}"#)
            })
            .collect();
        format!(r#","attributes":[{}]"#, items.join(","))
    }

    /// A cache directory holding one chromosome directory per `(chr, body)`
    /// pair, each with one shard of that body.
    fn cache_with_shards(prefix: &str, shards: &[(&str, &str)]) -> tempfile::TempDir {
        let dir = tempfile::Builder::new().prefix(prefix).tempdir().unwrap();
        for (chr, body) in shards {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            std::fs::write(chr_dir.join("1-1000000.json"), body).unwrap();
        }
        dir
    }

    /// `TSL` follows VEP's `m/tsl(\d+)/` with `$hash->{TSL} = $1 if $1`
    /// (OutputFactory.pm `BaseTranscriptVariationAllele_to_output_hash`, the
    /// same rule the pick order reads): `tsl1 (assigned to previous version 8)`
    /// is 1, `tslNA` with or without a suffix is none, and a zero is none.
    #[test]
    fn tsl_leading_digits_with_assigned_suffix() {
        let tsl = |value: &str| transcript_from(&minimal(&attributes(&[("TSL", value)]))).tsl;
        assert_eq!(tsl("tsl1 (assigned to previous version 8)"), Some(1));
        assert_eq!(tsl("tsl5"), Some(5));
        assert_eq!(tsl("tsl3 (assigned to previous version 2)"), Some(3));
        assert_eq!(tsl("tslNA"), None);
        assert_eq!(tsl("tslNA (assigned to previous version 4)"), None);
        assert_eq!(tsl("tsl0"), None);
        // The first digit-bearing match is the only one read, as the regex's is.
        assert_eq!(tsl("tsl0 (assigned to previous version tsl3)"), None);
        assert_eq!(tsl("tslNA then tsl4"), Some(4));
        assert_eq!(tsl_from_attribute_value("tsl2"), Some(2));
        assert_eq!(tsl_from_attribute_value(""), None);

        // A cache writing the level as a top-level number is read as such.
        assert_eq!(transcript_from(&minimal(r#","tsl":2"#)).tsl, Some(2));
        assert_eq!(transcript_from(&minimal("")).tsl, None);
    }

    /// `appris` is the attribute's value as the cache spells it: VEP rewrites
    /// `principal` to `P` and `alternative` to `A` only when it prints the
    /// column, and its pick order reads the spelled value's letter and digit.
    #[test]
    fn appris_derived_from_attribute() {
        let appris =
            |value: &str| transcript_from(&minimal(&attributes(&[("appris", value)]))).appris;
        assert_eq!(appris("principal1").as_deref(), Some("principal1"));
        assert_eq!(appris("alternative2").as_deref(), Some("alternative2"));
        assert_eq!(appris("principalM").as_deref(), Some("principalM"));
        assert_eq!(transcript_from(&minimal("")).appris, None);
        // A top-level key is read before the attribute.
        let both = minimal(&format!(
            r#","appris":"principal3"{}"#,
            attributes(&[("appris", "alternative1")])
        ));
        assert_eq!(transcript_from(&both).appris.as_deref(), Some("principal3"));
    }

    /// The `MANE_Select` attribute's value, VEP's `MANE_SELECT` column, is a
    /// RefSeq accession; the field is `Some` on the attribute's presence,
    /// which is what lists `MANE_Select` under `MANE`.
    #[test]
    fn mane_select_derived_from_attribute() {
        let tx = transcript_from(&minimal(&attributes(&[("MANE_Select", "NM_198691.3")])));
        assert_eq!(tx.mane_select.as_deref(), Some("NM_198691.3"));
        assert_eq!(tx.mane_plus_clinical, None);
        assert_eq!(transcript_from(&minimal("")).mane_select, None);
        let empty = transcript_from(&minimal(&attributes(&[("MANE_Select", "")])));
        assert_eq!(empty.mane_select.as_deref(), Some(""));
    }

    /// `MANE_Plus_Clinical` is read beside `MANE_Select` from its own attribute
    /// (OutputFactory.pm lists `MANE_Plus_Clinical` under `MANE` on its presence
    /// and prints the value as `MANE_PLUS_CLINICAL`); a transcript may carry
    /// both.
    #[test]
    fn mane_plus_clinical_derived_from_attribute() {
        let tx = transcript_from(&minimal(&attributes(&[
            ("MANE_Plus_Clinical", "NM_001005484.2"),
            ("MANE_Select", "NM_001005484.3"),
        ])));
        assert_eq!(tx.mane_plus_clinical.as_deref(), Some("NM_001005484.2"));
        assert_eq!(tx.mane_select.as_deref(), Some("NM_001005484.3"));
        let select_only = transcript_from(&minimal(&attributes(&[("MANE_Select", "NM_1.1")])));
        assert_eq!(select_only.mane_plus_clinical, None);
        // The first attribute of a code is the one read, as VEP's grep takes it.
        let twice = transcript_from(&minimal(&attributes(&[
            ("MANE_Plus_Clinical", "NM_2.1"),
            ("MANE_Plus_Clinical", "NM_3.1"),
        ])));
        assert_eq!(twice.mane_plus_clinical.as_deref(), Some("NM_2.1"));
    }

    /// GENCODE set membership is the presence of the `gencode_basic` and
    /// `gencode_primary` attributes, whose dumped values are the labels
    /// `GENCODE basic` and `GENCODE primary` (AnnotationType/Transcript.pm
    /// `filter_transcript` greps the codes alone); one set does not imply the
    /// other.
    #[test]
    fn gencode_basic_and_primary_by_presence() {
        let both = transcript_from(&minimal(&attributes(&[
            ("gencode_basic", "GENCODE basic"),
            ("gencode_primary", "GENCODE primary"),
        ])));
        assert!(both.gencode_basic && both.gencode_primary);
        let basic = transcript_from(&minimal(&attributes(&[("gencode_basic", "GENCODE basic")])));
        assert!(basic.gencode_basic && !basic.gencode_primary);
        let numeric = transcript_from(&minimal(&attributes(&[("gencode_primary", "1")])));
        assert!(!numeric.gencode_basic && numeric.gencode_primary);
        let neither = transcript_from(&minimal(&attributes(&[("TSL", "tsl1")])));
        assert!(!neither.gencode_basic && !neither.gencode_primary);
    }

    /// Every chromosome directory under `transcripts/` is a chromosome of the
    /// cache, patches, alternate loci, scaffolds and `LRG_` regions included, in
    /// the sorted order of VEP's `valid_chromosomes` (CacheDir.pm reads the
    /// cache directory's subdirectory names); each transcript is stamped with its
    /// directory's name.
    #[test]
    fn alt_contig_directory_loads() {
        let tx = |id: &str| {
            format!(
                r#"[{{"stable_id":"{id}","gene_stable_id":"ENSG_{id}","start":100,"end":900,"strand":1,"biotype":"protein_coding","source":"ensembl"}}]"#
            )
        };
        let dir = cache_with_shards(
            "vep_jc_alt_",
            &[
                ("HG1012_PATCH", &tx("ENST_PATCH")),
                ("21", &tx("ENST_21")),
                ("GL000008.2", &tx("ENST_GL")),
                ("LRG_1", &tx("ENST_LRG")),
                ("KI270728.1", &tx("ENST_KI")),
            ],
        );
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let mut chromosomes: Vec<&String> = shards.keys().collect();
        chromosomes.sort();
        assert_eq!(
            chromosomes,
            ["21", "GL000008.2", "HG1012_PATCH", "KI270728.1", "LRG_1"]
        );
        for (chr, id) in [("HG1012_PATCH", "ENST_PATCH"), ("LRG_1", "ENST_LRG")] {
            let txs = load_transcripts_for_chr(chr, &shards[chr]).unwrap();
            assert_eq!(txs.len(), 1);
            assert_eq!(&*txs[0].stable_id, id);
            assert_eq!(txs[0].chr, chr);
        }
    }

    /// `is_canonical` is the dumped flag (Database/Transcript.pm sets it to 1
    /// on the gene's canonical transcript and leaves it undefined otherwise),
    /// read as an integer or a boolean.
    #[test]
    fn is_canonical_from_int_or_bool() {
        assert!(transcript_from(&minimal(r#","is_canonical":1"#)).canonical);
        assert!(transcript_from(&minimal(r#","is_canonical":true"#)).canonical);
        assert!(!transcript_from(&minimal(r#","is_canonical":0"#)).canonical);
        assert!(!transcript_from(&minimal(r#","is_canonical":false"#)).canonical);
        assert!(!transcript_from(&minimal("")).canonical);
    }

    /// The gene identifiers VEP prefetches under `_gene_symbol`,
    /// `_gene_symbol_source` and `_gene_hgnc_id` arrive under the converter's
    /// unprefixed names; a transcript of an unnamed gene carries none.
    #[test]
    fn gene_symbol_loaded_from_renamed_key() {
        let tx = transcript_from(&minimal(
            r#","gene_symbol":"KRTAP10-1","gene_symbol_source":"HGNC","gene_hgnc_id":"HGNC:22966""#,
        ));
        assert_eq!(tx.gene_symbol.as_deref(), Some("KRTAP10-1"));
        assert_eq!(tx.gene_symbol_source.as_deref(), Some("HGNC"));
        assert_eq!(tx.hgnc_id.as_deref(), Some("HGNC:22966"));
        let unnamed = transcript_from(&minimal(""));
        assert!(unnamed.gene_symbol.is_none() && unnamed.gene_symbol_source.is_none());
        assert!(unnamed.hgnc_id.is_none());
    }

    /// `gene_stable_id` (VEP's `_gene_stable_id`) is required: a shard whose
    /// transcript lacks it fails to parse, and the loaded transcript carries
    /// it; the gene's `_gene` object is not in the cache and nothing reads it.
    #[test]
    fn gene_stable_id_required() {
        assert_eq!(
            &*transcript_from(&minimal("")).gene_stable_id,
            "ENSG00000215455"
        );
        let without = r#"{"stable_id":"ENST1","start":1,"end":9,"strand":1,"biotype":"lncRNA","source":"ensembl"}"#;
        let err = serde_json::from_str::<JsonTranscript>(without)
            .err()
            .unwrap();
        assert!(err.to_string().contains("gene_stable_id"), "{err}");
    }

    /// `biotype` is required and read verbatim.
    #[test]
    fn biotype_required() {
        assert_eq!(&*transcript_from(&minimal("")).biotype, "protein_coding");
        let without = r#"{"stable_id":"ENST1","gene_stable_id":"ENSG1","start":1,"end":9,"strand":1,"source":"ensembl"}"#;
        let err = serde_json::from_str::<JsonTranscript>(without)
            .err()
            .unwrap();
        assert!(err.to_string().contains("biotype"), "{err}");
    }

    /// A Storable-derived transcript carries no chromosome of its own (the
    /// dumper keeps the slice, the converter drops it), so the shard directory
    /// names it; a `chr` key is kept when present. Coordinates are read from
    /// numbers or strings, and a negative strand is the reverse strand.
    #[test]
    fn chr_falls_back_to_shard_directory() {
        let tx = transcript_from(&minimal(""));
        assert_eq!(tx.chr, "21");
        assert_eq!((tx.start, tx.end), (44538981, 44540195));
        assert_eq!(tx.strand, Strand::Reverse);
        let named = transcript_from(
            r#"{"stable_id":"ENST1","gene_stable_id":"ENSG1","chr":"X","start":"100","end":"900","strand":"1","biotype":"lncRNA","source":"ensembl"}"#,
        );
        assert_eq!(named.chr, "X");
        assert_eq!((named.start, named.end), (100, 900));
        assert_eq!(named.strand, Strand::Forward);
    }

    /// The transcript `version` arrives as a number or a string holding one.
    #[test]
    fn version_from_number_or_string() {
        assert_eq!(
            transcript_from(&minimal(r#","version":3"#)).version,
            Some(3)
        );
        assert_eq!(
            transcript_from(&minimal(r#","version":"3""#)).version,
            Some(3)
        );
        assert_eq!(transcript_from(&minimal("")).version, None);
        assert_eq!(transcript_from(&minimal(r#","version":"""#)).version, None);
    }

    /// `gene_version` (VEP's `_gene_version`, the `.version` suffix of the Gene
    /// column under `--gene_version`) arrives as a number or a string.
    #[test]
    fn gene_version_loaded() {
        assert_eq!(
            transcript_from(&minimal(r#","gene_version":4"#)).gene_version,
            Some(4)
        );
        assert_eq!(
            transcript_from(&minimal(r#","gene_version":"12""#)).gene_version,
            Some(12)
        );
        assert_eq!(transcript_from(&minimal("")).gene_version, None);
    }

    /// `uniprot_isoform` (VEP's `_uniprot_isoform`) is read as the cache writes
    /// it, the `-` VEP stores for a translation without one included; the
    /// column suppresses the sentinel.
    #[test]
    fn uniprot_isoform_loaded() {
        assert_eq!(
            transcript_from(&minimal(r#","uniprot_isoform":"P60331-2""#))
                .uniprot_isoform
                .as_deref(),
            Some("P60331-2")
        );
        assert_eq!(
            transcript_from(&minimal(r#","uniprot_isoform":"-""#))
                .uniprot_isoform
                .as_deref(),
            Some("-")
        );
        assert_eq!(transcript_from(&minimal("")).uniprot_isoform, None);
    }

    /// The VEFC intron list (VEP's `get_all_Introns`) is the transcript's intron
    /// list when present, even where it disagrees with the gaps between exons;
    /// without one the gaps are derived.
    #[test]
    fn vefc_introns_preferred_over_derived() {
        let exons = r#","exons":[{"start":100,"end":200},{"start":301,"end":400}]"#;
        let with_vefc = transcript_from(&minimal(&format!(
            r#"{exons},"variation_effect_feature_cache":{{"introns":[{{"start":"201","end":"299"}}]}}"#
        )));
        let introns: Vec<(u64, u64)> = with_vefc.introns.iter().map(|i| (i.start, i.end)).collect();
        assert_eq!(introns, vec![(201, 299)]);
        assert_eq!(with_vefc.vefc.unwrap().introns.len(), 1);

        let derived = transcript_from(&minimal(exons));
        let introns: Vec<(u64, u64)> = derived.introns.iter().map(|i| (i.start, i.end)).collect();
        assert_eq!(introns, vec![(201, 300)]);
    }

    /// `peptide` (the translation's sequence with its edits applied) loads on a
    /// translated transcript and is absent on one without a translation.
    #[test]
    fn peptide_loaded_when_translated() {
        let coding = transcript_from(&minimal(
            r#","translation":{"stable_id":"ENSP00000383226","start":46,"end":894},"variation_effect_feature_cache":{"peptide":"MAASTMSVC","translateable_seq":"ATGGCCGCG"}"#,
        ));
        assert_eq!(coding.vefc.unwrap().peptide.as_deref(), Some("MAASTMSVC"));
        let non_coding = transcript_from(&minimal(
            r#","variation_effect_feature_cache":{"translateable_seq":""}"#,
        ));
        assert_eq!(non_coding.vefc.unwrap().peptide, None);
    }

    /// The runtime-only keys VEP's `BaseTranscriptVariation` writes into a live
    /// transcript's cache, and any key the loader has no use for, are ignored
    /// rather than refused.
    #[test]
    fn unknown_vefc_keys_ignored() {
        let tx = transcript_from(&minimal(
            r#","description":"keratin associated protein 10-1","variation_effect_feature_cache":{"_exon_interval_tree":{},"_intron_numbers":"1-2","_has_frameshift_intron":0,"translation_md5":"4fb6e2a4930c8b5aa1ac86f88fbc545b","spliced_seq":"ATG","codon_table":"2","bam":{}}"#,
        ));
        assert_eq!(tx.vefc.unwrap().codon_table, 2);
    }

    /// A protein feature's `analysis` is its display label: a string, or the
    /// dumped `Bio::EnsEMBL::Analysis` object reduced to `_display_label`
    /// (Database/Transcript.pm keeps that key alone); an object without a label
    /// yields none.
    #[test]
    fn protein_feature_label_from_display_label() {
        let tx = transcript_from(&minimal(
            r#","variation_effect_feature_cache":{"protein_features":[
                {"start":"1","end":"282","hseqname":"AF-P60331-F1","analysis":"AFDB-ENSP mappings"},
                {"start":22,"end":282,"hseqname":"PF00080","analysis":{"_display_label":"Pfam"}},
                {"start":5,"end":9,"hseqname":"PTHR24416","analysis":{"display_label":"PANTHER"}},
                {"start":5,"end":9,"hseqname":"orphan","analysis":{}},
                {"start":5,"end":9,"hseqname":"blank","analysis":""}]}"#,
        ));
        let labels: Vec<Option<&str>> = tx
            .vefc
            .as_ref()
            .unwrap()
            .protein_features
            .iter()
            .map(|pf| pf.analysis.as_deref())
            .collect();
        assert_eq!(
            labels,
            vec![
                Some("AFDB-ENSP mappings"),
                Some("Pfam"),
                Some("PANTHER"),
                None,
                None
            ]
        );
    }

    /// A prediction matrix arrives as base64 of the Storable `matrix` blob: a
    /// gzip stream over `VEP` and 20 little-endian 16-bit cells per residue,
    /// each the category in its top two bits over the score times 1000
    /// (ProteinFunctionPredictionMatrix.pm `pack 'v'`); a cell of all ones is
    /// no prediction, and a tool's empty object is no matrix.
    #[test]
    fn storable_matrix_decodes() {
        use base64::Engine;
        use std::io::Write as _;
        use vep_core::prediction::AnalysisType;

        let mut raw = b"VEP".to_vec();
        raw.extend(std::iter::repeat_n(0xFFu8, 2 * 20 * 40));
        // Residue 2, amino acid C (column 1): SIFT category 1 (deleterious),
        // score 0.01.
        let cell = ((1u16) << 14 | 10).to_le_bytes();
        let offset = 3 + (20 + 1) * 2;
        raw[offset..offset + 2].copy_from_slice(&cell);
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&raw).unwrap();
        let blob = base64::engine::general_purpose::STANDARD.encode(gz.finish().unwrap());

        let tx = transcript_from(&minimal(&format!(
            r#","variation_effect_feature_cache":{{"protein_function_predictions":{{"sift":{{"analysis":"sift","translation_md5":"4fb6e2a4930c8b5aa1ac86f88fbc545b","predictions_data":"{blob}"}},"polyphen_humvar":{{}},"polyphen_humdiv":{{}}}}}}"#
        )));
        let predictions = tx.vefc.unwrap().protein_function_predictions.unwrap();
        let sift = predictions.sift.unwrap();
        assert_eq!(sift.analysis, "sift");
        assert_eq!(sift.peptide_length, 0);
        let hit = sift.lookup(2, b'C', AnalysisType::Sift).unwrap();
        assert_eq!((hit.label.as_str(), hit.score), ("deleterious", 0.01));
        assert!(sift.lookup(1, b'A', AnalysisType::Sift).is_none());
        assert!(predictions.polyphen_humvar.is_none());
        assert!(predictions.polyphen_humdiv.is_none());
    }

    /// A cache whose `transcripts` directory holds no shard is refused under
    /// the error contract, naming the serialiser `info.json` declares: a Sereal
    /// cache converts to no shard. A cache without an `info.json` is a Storable
    /// one.
    #[test]
    fn sereal_cache_without_shards_is_refused() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_sereal_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        std::fs::write(chr_dir.join("1-1000000.sereal"), b"=srl").unwrap();
        std::fs::write(
            dir.path().join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","serialiser_type":"sereal"}"#,
        )
        .unwrap();
        let root = dir.path().to_str().unwrap();
        let err = enumerate_transcript_shards(root).unwrap_err();
        assert_eq!(
            CliError::from_run(&err),
            format!(
                "ERROR: [cache] no transcript shards under {}/transcripts (serialiser_type sereal)",
                root
            )
        );

        std::fs::remove_file(dir.path().join("info.json")).unwrap();
        let err = enumerate_transcript_shards(root).unwrap_err();
        assert!(
            CliError::from_run(&err).ends_with("(serialiser_type storable)"),
            "{err}"
        );
    }

    /// The spliced cDNA of a coding transcript is assembled from the cached
    /// UTRs and translateable sequence as the transcript is read.
    #[test]
    fn spliced_seq_assembled_for_coding_transcript() {
        let tx = transcript_from(&minimal(
            r#","translation":{"stable_id":"ENSP1","start":4,"end":12},"exons":[{"start":44538981,"end":44538995}],
            "variation_effect_feature_cache":{"five_prime_utr":"GGA","translateable_seq":"ATGCCCTAA","three_prime_utr":"TTT",
            "mapper":{"start_phase":-1,"cdna_coding_start":4,"cdna_coding_end":12,"pairs":[{"from_start":1,"from_end":15,"to_start":44538981,"to_end":44538995,"ori":-1}]}}"#,
        ));
        assert_eq!(tx.spliced_seq(), Some("GGAATGCCCTAATTT"));
        assert_eq!(transcript_from(&minimal("")).spliced_seq(), None);
    }

    /// The variation record's ClinVar allele fields and the gnomAD genomes
    /// Amish frequency load under VEP's column names; an empty string is
    /// absent, and a record with none of the three carries no ClinVar group.
    #[test]
    fn variation_clinical_fields_loaded() {
        let json = r#"{"variation_name":"rs1","start":100,"allele_string":"A/G","strand":1,
            "clin_sig":"pathogenic","clin_sig_allele":"G:pathogenic;G:likely_pathogenic","clin_sig_ref_allele":"A",
            "clinical_impact":"Melanoma(Tier I - Strong);Melanoma(oncogenicity:Oncogenic)","gnomADg_AMI":"G:0.0123"}"#;
        let cv = convert_variation(serde_json::from_str(json).unwrap());
        let clinvar = cv.clinvar.as_deref().unwrap();
        assert_eq!(
            clinvar.clin_sig_allele.as_deref(),
            Some("G:pathogenic;G:likely_pathogenic")
        );
        assert_eq!(clinvar.clin_sig_ref_allele.as_deref(), Some("A"));
        assert_eq!(
            clinvar.clinical_impact.as_deref(),
            Some("Melanoma(Tier I - Strong);Melanoma(oncogenicity:Oncogenic)")
        );
        assert_eq!(
            cv.frequencies.get("gnomADg_AMI").map(String::as_str),
            Some("G:0.0123")
        );

        let empty = r#"{"variation_name":"rs2","start":100,"allele_string":"A/G","clin_sig_allele":"","clin_sig_ref_allele":"","clinical_impact":"","gnomADg_AMI":""}"#;
        let cv = convert_variation(serde_json::from_str(empty).unwrap());
        assert!(cv.clinvar.is_none());
        assert!(!cv.frequencies.contains_key("gnomADg_AMI"));
    }
}
