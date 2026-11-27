// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Load transcripts from the converted JSON cache directory.
//!
//! The JSON cache (written by `vep-cache-builder` or converted from an Ensembl VEP cache) uses field names that
//! differ from the Rust `Transcript` struct (e.g. `dbID` vs `db_id`,
//! `is_canonical` vs `canonical`). This module provides deserialization types
//! that map the JSON format and converts them to the canonical Rust types.
//!
//! Ensembl VEP source citations name modules of ensembl-vep (`Bio/EnsEMBL/VEP/...`);
//! `VariationEffect.pm` is `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, and
//! `BaseTranscriptVariation.pm`, `BaseVariationFeature.pm` and
//! `ProteinFunctionPredictionMatrix.pm` are
//! `Bio/EnsEMBL/Variation/...`, in ensembl-variation; `TranscriptMapper.pm` is
//! `Bio/EnsEMBL/TranscriptMapper.pm` in ensembl core.

use std::cell::Cell;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::Deserialize;
use tracing::info;

use vep_core::cache_info::CacheInfo;
use vep_core::cache_manifest::CacheManifest;
use vep_core::coordinate::Strand;
use vep_core::transcript::{
    tsl_level, Attribute, BamEditStatus, Exon, ExonCoordMapper, Intron, MapperPair,
    PredictionMatrix, ProteinFeature, ProteinFunctionPredictions, SeqEdit, SourceCache, Transcript,
    TranscriptMapper, TranscriptVEFC, Translation,
};

use crate::cache_resolve::manifest_error;
use crate::config::CliError;
use crate::transcript_index::TranscriptSet;

/// Visitors for the numeric fields a Storable-derived cache may store as a JSON number
/// or as a string holding one.
///
/// An `#[serde(untagged)]` enum would buffer each value and then try its
/// variants in turn, formatting an error message for every variant that does
/// not match, so a plain integer field cost an allocation and a formatted
/// error on the load path. Each visitor below accepts the same inputs with the
/// same results: an integer as itself, a negative integer as 0, a float
/// truncated, a string parsed (an empty string is "absent" where the field is
/// optional, and a negative integer in a string is 0, as the number is), and
/// for `u8` a bool as 0/1.
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
                    return Ok(None);
                }
                match v.parse::<$t>() {
                    Ok(n) => Ok(Some(n)),
                    Err(error) => match v.parse::<i64>() {
                        Ok(n) if n < 0 => Ok(Some(0)),
                        _ => Err(E::custom(error)),
                    },
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

/// Visitor for a genomic coordinate that can lie below 1, as the start of a
/// transcript stored across the origin of a circular sequence does
/// ([`place_across_origin`]): the inputs of [`U64OrStr`], a negative integer
/// kept as it is.
struct SignedOrStr;
impl serde::de::Visitor<'_> for SignedOrStr {
    type Value = Option<i64>;

    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("an integer, a string holding one, or an empty string")
    }

    fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
        Ok(Some(i64::try_from(v).unwrap_or(i64::MAX)))
    }

    fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
        Ok(Some(v))
    }

    fn visit_f64<E: serde::de::Error>(self, v: f64) -> Result<Self::Value, E> {
        Ok(Some(v as i64))
    }

    fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
        if v.is_empty() {
            Ok(None)
        } else {
            v.parse::<i64>().map(Some).map_err(E::custom)
        }
    }
}

fn deserialize_signed_string_or_number<'de, D>(deserializer: D) -> Result<i64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer
        .deserialize_any(SignedOrStr)?
        .ok_or_else(|| serde::de::Error::custom("cannot parse integer from empty string"))
}

/// A stored coordinate as a position on its sequence: one below 1 reads as 0.
fn position(stored: i64) -> u64 {
    u64::try_from(stored).unwrap_or(0)
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

/// Visitor for a transcript's top-level `tsl` key, a support level a cache may
/// write beside the `TSL` attribute: a non-negative integer or a string of
/// digits, carried as the digit run's text; an empty string is none.
struct LevelOrStr;
impl serde::de::Visitor<'_> for LevelOrStr {
    type Value = Option<String>;

    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("a non-negative integer, a string of digits, or an empty string")
    }

    fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
        Ok(Some(v.to_string()))
    }

    fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
        u64::try_from(v)
            .map(|level| Some(level.to_string()))
            .map_err(|_| E::invalid_value(serde::de::Unexpected::Signed(v), &self))
    }

    fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
        if v.is_empty() {
            Ok(None)
        } else if v.bytes().all(|b| b.is_ascii_digit()) {
            Ok(Some(v.to_string()))
        } else {
            Err(E::invalid_value(serde::de::Unexpected::Str(v), &self))
        }
    }
}

fn deserialize_optional_level<'de, D>(deserializer: D) -> Result<Option<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    deserializer.deserialize_option(OptionOf(LevelOrStr))
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
pub(crate) fn deserialize_arc_str<'de, D>(deserializer: D) -> Result<Arc<str>, D::Error>
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

pub(crate) fn deserialize_optional_arc_str<'de, D>(
    deserializer: D,
) -> Result<Option<Arc<str>>, D::Error>
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
    /// Absent on a transcript whose gene has no stable id, which Ensembl VEP
    /// keeps and annotates (`AnnotationType/Transcript.pm:220` and `:280`
    /// test the transcript's own id alone).
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_stable_id: Option<Arc<str>>,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    gene_version: Option<u32>,
    #[serde(default)]
    chr: Option<String>,
    #[serde(deserialize_with = "deserialize_signed_string_or_number")]
    start: i64,
    #[serde(deserialize_with = "deserialize_string_or_number")]
    end: u64,
    #[serde(deserialize_with = "deserialize_strand")]
    strand: Strand,
    #[serde(deserialize_with = "deserialize_arc_str")]
    biotype: Arc<str>,
    source: String,
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_symbol: Option<Arc<str>>,
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_symbol_source: Option<Arc<str>>,
    #[serde(default, deserialize_with = "deserialize_optional_arc_str")]
    gene_symbol_id: Option<Arc<str>>,
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
    #[serde(default, deserialize_with = "deserialize_optional_level")]
    tsl: Option<String>,
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
    source_cache: Option<String>,
    #[serde(default)]
    bam_edit_status: Option<String>,
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

thread_local! {
    /// Whether the shard parse running on this thread keeps each exon's stable
    /// id ([`Transcript::exon_ids`]): off, the shard's `stable_id` key is
    /// consumed and dropped, since no output or consequence reads an exon's id
    /// and parsing it is a cost of every load; on, for a run whose plugins print
    /// the ids (NearestExonJB). Set by [`parse_shard_transcripts`] for the
    /// duration of its parse, so the loads of one process read their own option.
    static KEEP_EXON_IDS: Cell<bool> = const { Cell::new(false) };
    /// Whether the shard parse running on this thread keeps each transcript's
    /// SIFT and PolyPhen matrices: off for a run without `--sift` or
    /// `--polyphen`, which never reads them and whose shards are mostly their
    /// bytes. Set by [`parse_shard_transcripts`] for the duration of its parse.
    static KEEP_PREDICTIONS: Cell<bool> = const { Cell::new(true) };
}

/// A transcript's `protein_function_predictions` key: read when the parse keeps
/// the matrices, else consumed and dropped.
fn deserialize_predictions<'de, D>(
    deserializer: D,
) -> Result<Option<JsonProteinFunctionPredictions>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    if KEEP_PREDICTIONS.get() {
        Option::deserialize(deserializer)
    } else {
        serde::de::IgnoredAny::deserialize(deserializer).map(|_| None)
    }
}

/// An exon's `stable_id` key: read as an id when the parse keeps them, else
/// consumed and dropped.
fn deserialize_exon_id<'de, D>(deserializer: D) -> Result<Option<Arc<str>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    if KEEP_EXON_IDS.get() {
        deserialize_optional_arc_str(deserializer)
    } else {
        serde::de::IgnoredAny::deserialize(deserializer).map(|_| None)
    }
}

/// An exon's span, rank and phases, and its stable id when the parse keeps
/// exon ids.
#[derive(Deserialize)]
struct JsonExon {
    #[serde(default, deserialize_with = "deserialize_exon_id")]
    stable_id: Option<Arc<str>>,
    #[serde(deserialize_with = "deserialize_signed_string_or_number")]
    start: i64,
    #[serde(deserialize_with = "deserialize_signed_string_or_number")]
    end: i64,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    rank: Option<u32>,
    #[serde(default, deserialize_with = "deserialize_i8_from_number")]
    phase: i8,
    #[serde(default, deserialize_with = "deserialize_i8_from_number")]
    end_phase: i8,
}

#[derive(Deserialize)]
struct JsonIntron {
    #[serde(deserialize_with = "deserialize_signed_string_or_number")]
    start: i64,
    #[serde(deserialize_with = "deserialize_signed_string_or_number")]
    end: i64,
    #[serde(default, deserialize_with = "deserialize_optional_u32")]
    rank: Option<u32>,
}

/// A transcript's translation. The converter writes defined values alone, so
/// a translation the cache holds without a stable id has no `stable_id` key and
/// loads with an empty id. Ensembl VEP reads that id in one place, the
/// reference name of an HGVSp string (`TranscriptVariationAllele.pm:1748` at
/// ensembl-variation release/116).
#[derive(Deserialize)]
struct JsonTranslation {
    #[serde(default)]
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
    spliced_seq: Option<String>,
    #[serde(default)]
    introns: Vec<JsonIntron>,
    #[serde(default)]
    sorted_exons: Vec<JsonExon>,
    #[serde(default)]
    mapper: Option<JsonMapper>,
    #[serde(default)]
    protein_features: Vec<JsonProteinFeature>,
    #[serde(default, deserialize_with = "deserialize_predictions")]
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
    /// The gzipped prediction matrix, base64 in the shard and decoded here
    /// straight from the shard's bytes; empty when the key is absent or its
    /// text is not base64.
    #[serde(default, deserialize_with = "deserialize_base64_bytes")]
    predictions_data: Vec<u8>,
}

/// Decode a base64 string field into its bytes without an intermediate
/// `String`. Text that is not base64 decodes to no bytes, which
/// [`convert_prediction_matrix`] reads as no matrix.
fn deserialize_base64_bytes<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use base64::Engine;

    struct Base64Visitor;
    impl serde::de::Visitor<'_> for Base64Visitor {
        type Value = Vec<u8>;

        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a base64 string")
        }

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
            Ok(base64::engine::general_purpose::STANDARD
                .decode(v)
                .unwrap_or_default())
        }
    }
    deserializer.deserialize_str(Base64Visitor)
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

/// The matrix a shard's record carries, its gzip stream kept as stored and
/// inflated on the first lookup ([`PredictionMatrix::from_compressed`]); a
/// record without a stream (a tool's empty object, or text that is not
/// base64) is no matrix.
fn convert_prediction_matrix(jpm: JsonPredictionMatrix) -> Option<PredictionMatrix> {
    if jpm.predictions_data.is_empty() {
        return None;
    }
    Some(PredictionMatrix::from_compressed(
        jpm.analysis,
        jpm.sub_analysis,
        jpm.peptide_length,
        jpm.predictions_data,
    ))
}

fn convert_protein_function_predictions(
    jp: JsonProteinFunctionPredictions,
) -> ProteinFunctionPredictions {
    ProteinFunctionPredictions {
        sift: jp.sift.and_then(convert_prediction_matrix),
        polyphen_humdiv: jp.polyphen_humdiv.and_then(convert_prediction_matrix),
        polyphen_humvar: jp.polyphen_humvar.and_then(convert_prediction_matrix),
    }
}

/// The exons' stable ids in rank order, `Some` when an exon carries an id (a
/// parse keeping exon ids), an exon without one holding the empty id; `None`
/// when no exon carries one: a parse dropping the ids, or a model whose exons
/// have none. The order of [`Transcript::exons`] once sorted by rank, so an
/// exon and its id share an index.
fn exon_ids_in_rank_order(exons: &[JsonExon]) -> Option<Arc<[Arc<str>]>> {
    if exons.iter().all(|e| e.stable_id.is_none()) {
        return None;
    }
    let mut ids: Vec<(u32, Arc<str>)> = exons
        .iter()
        .enumerate()
        .map(|(i, e)| {
            let id = e.stable_id.clone().unwrap_or_else(|| Arc::from(""));
            (e.rank.unwrap_or(i as u32 + 1), id)
        })
        .collect();
    ids.sort_by_key(|(rank, _)| *rank);
    Some(ids.into_iter().map(|(_, id)| id).collect())
}

fn convert_exon(je: JsonExon, rank_fallback: u32) -> Exon {
    Exon {
        start: position(je.start),
        end: position(je.end),
        rank: je.rank.unwrap_or(rank_fallback),
        phase: je.phase,
        end_phase: je.end_phase,
    }
}

fn convert_intron(ji: JsonIntron, rank_fallback: u32) -> Intron {
    Intron {
        start: position(ji.start),
        end: position(ji.end),
        rank: ji.rank.unwrap_or(rank_fallback),
    }
}

/// Transcript flags are the codes of every `cds_` attribute, whatever its value:
/// VEP prints them that way (OutputFactory.pm) and its consequence predicates
/// test the attribute's presence, not its value (VariationEffect.pm
/// `_overlaps_start_codon`, `_overlaps_stop_codon`).
fn derive_flags_from_attributes(attributes: &[Attribute]) -> Vec<String> {
    let mut flags: Vec<String> = Vec::new();
    for attr in attributes {
        let code: &str = &attr.code;
        if code.starts_with("cds_") && !flags.iter().any(|flag| flag == code) {
            flags.push(code.to_string());
        }
    }
    flags
}

/// The attribute codes of one shard, each allocated once: a cache repeats a
/// few dozen codes (`TSL`, `gencode_basic`, `appris`) on most of its
/// transcripts, so every attribute of the shard carrying a code shares the
/// one `Arc<str>` its first carrier made.
#[derive(Default)]
struct AttributeCodes(FxHashSet<Arc<str>>);

impl AttributeCodes {
    fn intern(&mut self, code: String) -> Arc<str> {
        if let Some(shared) = self.0.get(code.as_str()) {
            return Arc::clone(shared);
        }
        let shared: Arc<str> = Arc::from(code);
        self.0.insert(Arc::clone(&shared));
        shared
    }
}

/// The value of the first attribute carrying `code`: VEP reads MANE, TSL and
/// APPRIS from the attribute list with a `grep` on the code that keeps its
/// first match (OutputFactory.pm `BaseTranscriptVariationAllele_to_output_hash`),
/// so a later duplicate is never read.
fn attribute_value<'a>(attributes: &'a [Attribute], code: &str) -> Option<&'a str> {
    attributes
        .iter()
        .find(|a| &*a.code == code)
        .map(|a| a.value.as_str())
}

/// Whether any attribute carries `code`, whatever its value.
fn has_attribute(attributes: &[Attribute], code: &str) -> bool {
    attributes.iter().any(|a| &*a.code == code)
}

/// The `source_cache` key as the cache spells it (`Ensembl` or `RefSeq`,
/// `Pipeline/DumpVEP/MergeVEP.pm:103, 110`); any other value is refused, naming
/// the transcript, since the merged-cache filter and the `SOURCE` column would
/// read it as neither.
fn source_cache_of(jt: &JsonTranscript) -> Result<Option<SourceCache>, String> {
    match jt.source_cache.as_deref() {
        None => Ok(None),
        Some("Ensembl") => Ok(Some(SourceCache::Ensembl)),
        Some("RefSeq") => Ok(Some(SourceCache::RefSeq)),
        Some(other) => Err(format!(
            "transcript {}: source_cache `{other}` is neither Ensembl nor RefSeq",
            jt.stable_id
        )),
    }
}

/// The `bam_edit_status` key as `apply_edits` writes it (`ok` or `failed`,
/// `AnnotationType/Transcript.pm:601, 621`); any other value is refused, naming
/// the transcript.
fn bam_edit_status_of(jt: &JsonTranscript) -> Result<Option<BamEditStatus>, String> {
    match jt.bam_edit_status.as_deref() {
        None => Ok(None),
        Some("ok") => Ok(Some(BamEditStatus::Ok)),
        Some("failed") => Ok(Some(BamEditStatus::Failed)),
        Some(other) => Err(format!(
            "transcript {}: bam_edit_status `{other}` is neither ok nor failed",
            jt.stable_id
        )),
    }
}

/// Whether `set` admits a parsed record, read from the fields the set tests
/// before the transcript is built: the stable id, the chromosome the record
/// names or the shard's, `source_cache` (a value that is neither `Ensembl` nor
/// `RefSeq` is refused here, as the conversion refuses it) and the two GENCODE
/// attributes.
fn set_admits(
    set: &TranscriptSet,
    jt: &JsonTranscript,
    chr_override: &str,
) -> Result<bool, String> {
    let source_cache = source_cache_of(jt)?;
    let has = |code: &str| jt.attributes.iter().any(|attribute| attribute.code == code);
    Ok(set.accepts_model(
        &jt.stable_id,
        jt.chr.as_deref().unwrap_or(chr_override),
        source_cache.as_ref(),
        has("gencode_basic"),
        has("gencode_primary"),
    ))
}

/// Reads a transcript the cache stores across the origin of a circular
/// sequence into coordinates that run on past the sequence's last base, and
/// returns the sequence's length; `None`, and the record untouched, for every
/// other transcript.
///
/// Ensembl's dump stores such a transcript with a start below 1, its
/// database start less the sequence's length (ensembl core
/// `DBSQL/TranscriptAdaptor.pm:1941-1947`), and its exons the same way,
/// while the mapper pairs keep the database coordinates, the pair across the
/// origin with its start above its end (`TranscriptMapper.pm:160-161`). The
/// first base of the transcript among the pairs, less the stored start, is
/// therefore the sequence's length. Every coordinate at or below the stored
/// end lies among the sequence's first bases and is moved one length on, and
/// the pair across the origin has its end moved.
fn place_across_origin(jt: &mut JsonTranscript) -> Result<Option<u64>, String> {
    if jt.start >= 1 {
        return Ok(None);
    }
    let stored_end = jt.end;
    let first_base = jt
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.mapper.as_ref())
        .and_then(|mapper| {
            mapper
                .pairs
                .iter()
                .map(|pair| pair.to_start)
                .filter(|&start| start > stored_end)
                .min()
        })
        .ok_or_else(|| {
            format!(
                "transcript {}: its start {} lies below base 1 and no mapper pair places \
                 it across the origin of a circular sequence",
                jt.stable_id, jt.start
            )
        })?;
    // The stored start is the first base less the length, and is not above 0.
    let length = first_base.saturating_add(jt.start.unsigned_abs());
    let low_side = i64::try_from(stored_end).unwrap_or(i64::MAX);
    let turn = i64::try_from(length).unwrap_or(i64::MAX);
    let moved = |stored: i64| {
        if stored > low_side {
            stored
        } else {
            stored.saturating_add(turn)
        }
    };
    jt.start = moved(jt.start);
    jt.end = stored_end.saturating_add(length);
    for exon in &mut jt.exons {
        exon.start = moved(exon.start);
        exon.end = moved(exon.end);
    }
    if let Some(vefc) = jt.vefc.as_mut() {
        for exon in &mut vefc.sorted_exons {
            exon.start = moved(exon.start);
            exon.end = moved(exon.end);
        }
        for intron in &mut vefc.introns {
            intron.start = moved(intron.start);
            intron.end = moved(intron.end);
        }
        for pair in vefc.mapper.iter_mut().flat_map(|mapper| &mut mapper.pairs) {
            if pair.to_start <= stored_end {
                pair.to_start = pair.to_start.saturating_add(length);
                pair.to_end = pair.to_end.saturating_add(length);
            } else if pair.to_end <= stored_end
                && pair.from_start <= pair.from_end
                && length
                    .checked_sub(pair.to_start)
                    .and_then(|high| high.checked_add(pair.to_end))
                    == Some(pair.from_end - pair.from_start)
            {
                // The pair across the origin: its cDNA bases are those from its
                // start to the last base and from base 1 to its end.
                pair.to_end = pair.to_end.saturating_add(length);
            }
        }
    }
    Ok(Some(length))
}

/// The loaded transcript, or the reason a record's value cannot be represented.
/// `codes` shares one allocation among the attributes of a shard carrying the
/// same code.
fn convert_transcript(
    mut jt: JsonTranscript,
    chr_override: &str,
    codes: &mut AttributeCodes,
) -> Result<Transcript, String> {
    let circular_length = place_across_origin(&mut jt)?;
    let source_cache = source_cache_of(&jt)?;
    let bam_edit_status = bam_edit_status_of(&jt)?;
    let chr = jt.chr.unwrap_or_else(|| chr_override.to_string());

    let exon_ids = exon_ids_in_rank_order(&jt.exons);
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
        seq: t.seq,
    });

    let attributes: Vec<Attribute> = jt
        .attributes
        .into_iter()
        .map(|a| Attribute {
            code: codes.intern(a.code),
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
    // attributes; a top-level key, when a cache writes one, is read in the
    // attribute's place, and its one digit `0` is none, as the attribute's is
    // (`$1 if $1`, `OutputFactory.pm:1620` at release/116.3).
    let tsl = match jt.tsl {
        Some(level) => (level != "0").then_some(level),
        None => attribute_value(&attributes, "TSL")
            .and_then(tsl_level)
            .map(str::to_string),
    };
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

        let codon_table = v.codon_table.unwrap_or(1);
        let peptide = match (v.peptide, v.translateable_seq.as_deref()) {
            (Some(peptide), Some(coding)) => {
                Some(peptide_under_own_code(peptide, coding, codon_table))
            }
            (peptide, _) => peptide,
        };
        let mut vefc = TranscriptVEFC::default();
        vefc.codon_table = codon_table;
        vefc.five_prime_utr = v.five_prime_utr;
        vefc.three_prime_utr = v.three_prime_utr;
        vefc.translateable_seq = v.translateable_seq;
        vefc.peptide = peptide;
        vefc.spliced_seq = v.spliced_seq;
        vefc.introns = vefc_introns;
        vefc.sorted_exons = sorted_exons;
        vefc.mapper = mapper;
        vefc.protein_features = protein_features;
        vefc.protein_function_predictions = protein_function_predictions;
        vefc.seq_edits = seq_edits;

        (Some(vefc), tx_introns)
    } else {
        (None, derive_introns(&exons))
    };

    Ok({
        let mut transcript = Transcript::default();
        transcript.stable_id = jt.stable_id;
        transcript.version = jt.version;
        transcript.db_id = jt.db_id;
        transcript.gene_stable_id = jt.gene_stable_id.unwrap_or_default();
        transcript.gene_version = jt.gene_version;
        transcript.chr = chr;
        transcript.start = position(jt.start);
        transcript.end = jt.end;
        transcript.strand = jt.strand;
        transcript.biotype = jt.biotype;
        transcript.source = jt.source;
        transcript.gene_symbol = jt.gene_symbol;
        transcript.gene_symbol_source = jt.gene_symbol_source;
        transcript.gene_symbol_id = jt.gene_symbol_id;
        transcript.hgnc_id = jt.hgnc_id;
        transcript.gene_phenotype = jt.gene_phenotype;
        transcript.canonical = jt.canonical;
        transcript.mane_select = mane_select;
        transcript.mane_plus_clinical = mane_plus_clinical;
        transcript.tsl = tsl;
        transcript.appris = appris;
        transcript.ccds = jt.ccds;
        transcript.protein_id = jt.protein_id;
        transcript.refseq = jt.refseq;
        transcript.swissprot = jt.swissprot;
        transcript.trembl = jt.trembl;
        transcript.uniparc = jt.uniparc;
        transcript.uniprot_isoform = jt.uniprot_isoform;
        transcript.exons = exons;
        transcript.exon_ids = exon_ids;
        transcript.introns = introns;
        transcript.cdna_coding_start = jt.cdna_coding_start;
        transcript.cdna_coding_end = jt.cdna_coding_end;
        transcript.coding_region_start = coding_region_start;
        transcript.coding_region_end = coding_region_end;
        transcript.translation_start = coding_region_start;
        transcript.translation_end = coding_region_end;
        transcript.translation = translation;
        transcript.cdna_sequence = None;
        transcript.protein_sequence = None;
        transcript.flags = flags.into();
        transcript.gencode_basic = gencode_basic;
        transcript.gencode_primary = gencode_primary;
        transcript.attributes = attributes;
        transcript.source_cache = source_cache;
        transcript.bam_edit_status = bam_edit_status;
        transcript.circular_length = circular_length;
        transcript.region_files = None;
        transcript.vefc = vefc;
        transcript.derived = Default::default();
        transcript.id_rank = 0;
        transcript
    })
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

/// Compute coding region start/end from the mapper: the genomic positions of
/// the first and the last coding base a pair holds, in genomic order.
///
/// A coding bound can lie in bases the model carries and the genome lacks, a
/// block an `_rna_edit` inserts and no pair holds; the bound is then the
/// nearest coding base the genome has. Ensembl VEP takes both bounds from the
/// translation's exons (ensembl core `Transcript.pm:1087-1154`), which a
/// block leaves alone, and numbers coding positions from the cDNA coding
/// start (ensembl core `TranscriptMapper.pm:447-448`), so a record on such a
/// model keeps its coding columns. A coding sequence with no base on the
/// genome has no bounds.
fn compute_coding_region(mapper: &JsonMapper) -> (Option<u64>, Option<u64>) {
    let cds_start = match mapper.cdna_coding_start {
        Some(v) if v > 0 => v,
        _ => return (None, None),
    };
    let cds_end = match mapper.cdna_coding_end {
        Some(v) if v > 0 => v,
        _ => return (None, None),
    };
    // A coding end before the coding start, both on the cDNA, has the bounds
    // of the span between the two.
    let (cds_start, cds_end) = (cds_start.min(cds_end), cds_start.max(cds_end));

    let genomic = |pair: &JsonMapperPair, cdna: u64| {
        let offset = cdna - pair.from_start;
        if pair.ori >= 0 {
            pair.to_start + offset
        } else {
            pair.to_end - offset
        }
    };
    let coding = || {
        mapper.pairs.iter().filter(|pair| {
            pair.from_start <= pair.from_end
                && pair.from_end >= cds_start
                && pair.from_start <= cds_end
        })
    };
    let first = coding().min_by_key(|pair| pair.from_start);
    let last = coding().max_by_key(|pair| pair.from_end);
    let (Some(first), Some(last)) = (first, last) else {
        return (None, None);
    };
    let genomic_start = genomic(first, cds_start.max(first.from_start));
    let genomic_end = genomic(last, cds_end.min(last.from_end));

    // For reverse strand, start > end in genomic coords, so swap
    if genomic_start > genomic_end {
        (Some(genomic_end), Some(genomic_start))
    } else {
        (Some(genomic_start), Some(genomic_end))
    }
}

/// Enumerate the per-chromosome JSON shard paths in a cache directory.
///
/// Directory structure: `{json_cache_dir}/transcripts/{chr}/{start}-{end}.json`.
/// This is a `read_dir` walk only: no JSON is parsed and no transcript is
/// constructed, so it is cheap enough to run eagerly at startup even when the
/// run will only ever query one chromosome.
///
/// The returned map's key set is the complete set of chromosomes the cache
/// covers. That matters beyond convenience: several Ensembl VEP-parity predicates in
/// `runner.rs` / `annotator.rs` branch on whether a chromosome is *present* in
/// the cache (the `intergenic_variant` fallback, and the inter-chromosomal BND
/// mate-annotation gate), mirroring Ensembl VEP's own cache-loading behavior. A
/// lazy loader must therefore be able to answer "is this chromosome in the
/// cache?" without parsing it, which this function provides.
///
/// A cache directory without `transcripts/` is an error: an empty map would
/// annotate every variant as intergenic and exit 0. A `transcripts`
/// directory under which no chromosome holds a shard is refused under the error
/// contract, naming the cache's serialiser: a Sereal cache converts to no
/// shard at all, and a copy that landed no shard reads the same.
///
/// The key set is Ensembl VEP's `valid_chromosomes`: every subdirectory whose
/// name does not start with a dot (`CacheDir.pm:431`, `grep {!/^\./ && -d}`).
/// A chromosome directory holding no `.json` shard is a key with an empty shard
/// list, so a record on it is annotated `intergenic_variant` rather than
/// skipped as a chromosome the cache lacks; a dot-prefixed directory is no
/// chromosome, whatever it holds.
pub fn enumerate_transcript_shards(
    json_cache_dir: &str,
) -> anyhow::Result<HashMap<String, Vec<PathBuf>>> {
    let mut by_chr: HashMap<String, Vec<PathBuf>> = HashMap::new();
    let transcripts_dir = Path::new(json_cache_dir).join("transcripts");

    if !transcripts_dir.is_dir() {
        anyhow::bail!(
            "no cache found at {json_cache_dir}: no transcripts directory ({})",
            transcripts_dir.display()
        );
    }

    for chr_entry in std::fs::read_dir(&transcripts_dir)
        .with_context(|| format!("Failed to read {}", transcripts_dir.display()))?
    {
        let chr_entry = chr_entry?;
        // `path().is_dir()` follows a symlinked chromosome directory, which
        // `DirEntry::file_type` reports as a link.
        if !chr_entry.path().is_dir() {
            continue;
        }
        let chr = chr_entry.file_name().to_string_lossy().to_string();
        if chr.starts_with('.') {
            continue;
        }

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
        by_chr.insert(chr, shards);
    }

    if by_chr.values().all(Vec::is_empty) {
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
pub(crate) fn is_shard_path(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".json") || name.ends_with(".json.gz")
}

/// The numeric start of a shard's range from its file name; a name that does
/// not start with a number sorts last, in name order.
pub(crate) fn shard_start(path: &Path) -> (u64, String) {
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

/// The first and last position of the region a shard's file name states
/// (`<start>-<end>.json`, or the same under `.json.gz`); `None` for a name of
/// another shape.
pub(crate) fn shard_span(path: &Path) -> Option<(u64, u64)> {
    let name = path.file_name()?.to_str()?;
    let range = name
        .strip_suffix(".json.gz")
        .or_else(|| name.strip_suffix(".json"))?;
    let (start, end) = range.split_once('-')?;
    let (start, end) = (start.parse().ok()?, end.parse().ok()?);
    (start <= end).then_some((start, end))
}

/// Read a shard's bytes, checking them against the cache's manifest before a
/// `.json.gz` shard is inflated.
pub(crate) fn read_shard(path: &Path, manifest: &CacheManifest) -> anyhow::Result<Vec<u8>> {
    let raw = std::fs::read(path).with_context(|| format!("Failed to read {}", path.display()))?;
    manifest.verify_bytes(path, &raw).map_err(manifest_error)?;
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
/// vep-cache-builder and `tools/vep-storable-to-json` each write a bare array;
/// a shard that wraps the transcripts in a map of arrays
/// (`{"<key>": [transcript, null, ...]}`) is read too, its nulls dropped. The
/// first non-whitespace byte selects the shape, so a shard is parsed once and a
/// malformed shard reports the error from the parse that matched its shape. Map
/// keys are visited in sorted order so two loads of the same shard yield the
/// same transcript order.
///
/// Under `keep_exon_ids` each exon's `stable_id` is read
/// ([`Transcript::exon_ids`]), else consumed and dropped; the option holds for
/// this parse alone ([`KEEP_EXON_IDS`]), a panic inside the parse included.
/// Under `keep_predictions` the SIFT and PolyPhen matrices are read, else
/// consumed and dropped ([`KEEP_PREDICTIONS`]).
fn parse_shard_transcripts(
    data: &[u8],
    keep_exon_ids: bool,
    keep_predictions: bool,
) -> serde_json::Result<Vec<JsonTranscript>> {
    struct ResetOnExit;
    impl Drop for ResetOnExit {
        fn drop(&mut self) {
            KEEP_EXON_IDS.set(false);
            KEEP_PREDICTIONS.set(true);
        }
    }
    KEEP_EXON_IDS.set(keep_exon_ids);
    KEEP_PREDICTIONS.set(keep_predictions);
    let _reset = ResetOnExit;
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

/// The identity a loaded transcript is de-duplicated on: its source cache and
/// dbID, or without a dbID its stable id and version.
type TranscriptKey = (Option<SourceCache>, Option<u64>, Arc<str>, Option<u32>);

/// How many records a load parsed and how many of them carried the attribute
/// the run's GENCODE set tests ([`TranscriptSet::gencode_attribute`]); both
/// stay zero when the set names no GENCODE flag.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct GencodeTally {
    pub records: usize,
    pub carriers: usize,
}

impl GencodeTally {
    /// Adds another load's counts.
    pub fn add(&mut self, other: Self) {
        self.records += other.records;
        self.carriers += other.carriers;
    }

    /// The warning a load owes when `attribute` sat on none of the records it
    /// parsed: Ensembl VEP's `filter_transcript` drops each one
    /// (`AnnotationType/Transcript.pm:223-230`) and every record annotates as
    /// intergenic, which a cache of a species outside GENCODE does on every
    /// run under the flag. `None` when a record carried it or none was parsed.
    pub fn warning(&self, attribute: &str) -> Option<String> {
        (self.records > 0 && self.carriers == 0).then(|| {
            format!(
                "no transcript of the cache carries the {attribute} attribute; \
                 every transcript is filtered"
            )
        })
    }
}

/// The transcripts of one shard's bytes held in memory, each stamped with
/// `chr`: the records of a gene-model file read beside the cache
/// ([`crate::gene_model`]), which take the conversion and the checks a cache
/// shard's records take. No transcript set reads them, as Ensembl VEP's file
/// sources never call `filter_transcript` (`AnnotationType/Transcript.pm:215`
/// is the cache's and the database's alone), and none is ranked: the caller
/// ranks the chromosome once every source's transcripts are in
/// ([`assign_id_ranks`]).
pub(crate) fn transcripts_of_shard_bytes(
    chr: &str,
    data: &[u8],
    keep_exon_ids: bool,
) -> Result<Vec<Transcript>, String> {
    let json_transcripts =
        parse_shard_transcripts(data, keep_exon_ids, true).map_err(|error| error.to_string())?;
    let mut codes = AttributeCodes::default();
    json_transcripts
        .into_iter()
        .map(|jt| {
            let transcript = convert_transcript(jt, chr, &mut codes)?;
            check_mapper_bounds(&transcript)?;
            check_codon_table(&transcript)?;
            Ok(transcript)
        })
        .collect()
}

/// One chromosome's transcripts as a load returns them, with the load's
/// GENCODE tally.
#[derive(Debug)]
pub(crate) struct LoadedTranscripts {
    pub transcripts: Vec<Transcript>,
    pub gencode: GencodeTally,
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
/// transcript by `convert_transcript`. Every shard's bytes are checked against
/// `manifest`, the cache's checksum manifest, as they are read.
pub fn load_transcripts_for_chr(
    chr: &str,
    shards: &[PathBuf],
    manifest: &CacheManifest,
) -> anyhow::Result<Vec<Transcript>> {
    Ok(load_transcripts_in_set(chr, shards, manifest, None)?.transcripts)
}

/// [`load_transcripts_for_chr`] keeping the models of `set` alone, each dropped
/// record costing its parse and nothing else: the set reads the parsed record
/// before the transcript, its exons and its sequences are built, as Ensembl VEP
/// filters each region's records before it touches them
/// (`AnnotationSource/Cache/Transcript.pm:239-241`). A dropped record is never
/// checked, so a mapper or codon table the loader refuses in a kept record
/// passes unseen in a dropped one, as it does there; a `source_cache` value
/// the cache cannot hold is refused whichever way the set reads. Under a
/// GENCODE flag the records parsed and the records carrying the set's
/// attribute are counted ([`GencodeTally`]); the caller reads the count once
/// every load of the run is in, since a chromosome without a carrier says
/// nothing about the cache.
pub(crate) fn load_transcripts_in_set(
    chr: &str,
    shards: &[PathBuf],
    manifest: &CacheManifest,
    set: Option<&TranscriptSet>,
) -> anyhow::Result<LoadedTranscripts> {
    load_transcripts_in_set_with(chr, shards, manifest, set, false, true)
}

/// [`load_transcripts_in_set`] keeping each exon's stable id on its transcript
/// under `keep_exon_ids` ([`Transcript::exon_ids`]), for a run whose plugins
/// print the ids, and the SIFT and PolyPhen matrices under `keep_predictions`.
pub(crate) fn load_transcripts_in_set_with(
    chr: &str,
    shards: &[PathBuf],
    manifest: &CacheManifest,
    set: Option<&TranscriptSet>,
    keep_exon_ids: bool,
    keep_predictions: bool,
) -> anyhow::Result<LoadedTranscripts> {
    let options = ShardOptions {
        set,
        keep_exon_ids,
        keep_predictions,
    };
    let parsed = parse_shards(chr, shards, manifest, options);
    merge_shards(parsed, shards)
}

/// How [`parse_shards`] reads each shard: the run's transcript set and the
/// keys [`load_transcripts_in_set_with`] keeps.
#[derive(Clone, Copy)]
pub(crate) struct ShardOptions<'a> {
    pub(crate) set: Option<&'a TranscriptSet>,
    pub(crate) keep_exon_ids: bool,
    pub(crate) keep_predictions: bool,
}

/// One shard's transcripts as [`parse_shards`] reads them, before the
/// chromosome's copies are merged ([`merge_shards`]), with the lowest start of
/// any record the shard holds, the set's filter aside.
pub(crate) struct ParsedShard {
    transcripts: Vec<Transcript>,
    tally: GencodeTally,
    pub(crate) lowest_start: Option<i64>,
}

/// Parses `shards` in parallel on the global pool, each checked against
/// `manifest`, into the transcripts the run's set admits.
pub(crate) fn parse_shards(
    chr: &str,
    shards: &[PathBuf],
    manifest: &CacheManifest,
    options: ShardOptions<'_>,
) -> Vec<anyhow::Result<ParsedShard>> {
    let ShardOptions {
        set,
        keep_exon_ids,
        keep_predictions,
    } = options;
    let gencode_attribute = set.and_then(TranscriptSet::gencode_attribute);
    shards
        .par_iter()
        .map(|path| -> anyhow::Result<ParsedShard> {
            let data = read_shard(path, manifest)?;
            let json_transcripts = parse_shard_transcripts(&data, keep_exon_ids, keep_predictions)
                .with_context(|| format!("Failed to parse {}", path.display()))?;
            let tally = match gencode_attribute {
                Some(code) => GencodeTally {
                    records: json_transcripts.len(),
                    carriers: json_transcripts
                        .iter()
                        .filter(|jt| jt.attributes.iter().any(|attribute| attribute.code == code))
                        .count(),
                },
                None => GencodeTally::default(),
            };
            let lowest_start = json_transcripts.iter().map(|jt| jt.start).min();
            let mut codes = AttributeCodes::default();
            let mut transcripts: Vec<Transcript> = json_transcripts
                .into_iter()
                .filter_map(|jt| match set.map(|set| set_admits(set, &jt, chr)) {
                    Some(Ok(false)) => None,
                    Some(Err(message)) => Some(Err(message)),
                    Some(Ok(true)) | None => Some(convert_transcript(jt, chr, &mut codes)),
                })
                .collect::<Result<_, String>>()
                .map_err(|message| CliError::new("cache", message))?;
            for transcript in &transcripts {
                check_mapper_bounds(transcript)
                    .and_then(|()| check_codon_table(transcript))
                    .map_err(|message| CliError::new("cache", message))?;
            }
            copy_hgnc_ids_to_refseq_models(&mut transcripts);
            Ok(ParsedShard {
                transcripts,
                tally,
                lowest_start,
            })
        })
        .collect()
}

/// The chromosome's transcripts from the parses of `shards`, in the shards'
/// order: each transcript once, held to the region files it was read from, and
/// ranked.
pub(crate) fn merge_shards(
    parsed: Vec<anyhow::Result<ParsedShard>>,
    shards: &[PathBuf],
) -> anyhow::Result<LoadedTranscripts> {
    // A cache converted from VEP's region files stores a gene's transcripts in
    // every 1 Mb shard the gene overlaps, so the same transcript arrives from two
    // shards; VEP de-duplicates them by dbID on load, and so does this loader,
    // keeping the first copy. Caches without dbIDs are keyed on stable id and
    // version. A merged cache holds two dbID spaces, the Ensembl core's and the
    // RefSeq otherfeatures', so the key carries the source: a RefSeq model whose
    // dbID an Ensembl model shares elsewhere on the chromosome is a different
    // transcript, not a copy.
    let mut out: Vec<Transcript> = Vec::new();
    let mut gencode = GencodeTally::default();
    // Each kept transcript's place in `out`, and beside `out` the span of
    // every shard the transcript was read from.
    let mut seen: HashMap<TranscriptKey, usize> = HashMap::new();
    let mut read_from: Vec<Vec<Option<(u64, u64)>>> = Vec::new();
    for (item, shard) in parsed.into_iter().zip(shards) {
        let ParsedShard {
            transcripts, tally, ..
        } = item?;
        gencode.add(tally);
        let span = shard_span(shard);
        for t in transcripts {
            let key = match t.db_id {
                Some(id) => (t.source_cache.clone(), Some(id), Arc::from(""), None),
                None => (t.source_cache.clone(), None, t.stable_id.clone(), t.version),
            };
            match seen.entry(key) {
                std::collections::hash_map::Entry::Occupied(place) => {
                    read_from[*place.get()].push(span);
                }
                std::collections::hash_map::Entry::Vacant(place) => {
                    place.insert(out.len());
                    read_from.push(vec![span]);
                    out.push(t);
                }
            }
        }
    }
    hold_transcripts_to_their_region_files(&mut out, read_from);
    assign_id_ranks(&mut out);
    Ok(LoadedTranscripts {
        transcripts: out,
        gencode,
    })
}

/// Holds to the region files it was read from ([`Transcript::region_files`])
/// every transcript whose stable id another transcript of the chromosome
/// shares, and every transcript those files hold short of its span
/// (`filed_short_of_its_span`). A run of Ensembl VEP reads the files of the
/// regions a record reaches (`AnnotationSource.pm:190`), so a transcript is a
/// transcript of the regions whose files hold it.
///
/// Ensembl VEP's dump keeps one transcript a stable id in each region file of
/// a RefSeq or merged cache, the `ensembl`-sourced one, else the lowest `dbID`
/// (`Pipeline/DumpVEP/Dumper/Core.pm:109`,
/// `AnnotationType/Transcript.pm:313-343`). A transcript that lost its id to
/// another in one region's file and stands in a neighbour's is therefore no
/// transcript of the first region, though its span lies there.
///
/// `read_from` gives, beside each transcript, the span of every shard that
/// held it; a transcript read from a shard whose name states no span is left
/// to be found by its coordinates.
fn hold_transcripts_to_their_region_files(
    transcripts: &mut [Transcript],
    read_from: Vec<Vec<Option<(u64, u64)>>>,
) {
    let mut holders: HashMap<Arc<str>, u32> = HashMap::new();
    for transcript in transcripts.iter() {
        *holders.entry(transcript.stable_id.clone()).or_default() += 1;
    }
    for (transcript, spans) in transcripts.iter_mut().zip(read_from) {
        let Some(spans) = spans.into_iter().collect::<Option<Vec<_>>>() else {
            continue;
        };
        if holders[&transcript.stable_id] > 1 || filed_short_of_its_span(transcript, &spans) {
            transcript.region_files = Some(Arc::from(spans));
        }
    }
}

/// Whether the region files of `spans`, the ones `transcript` was read from,
/// leave a base of its span uncovered, on a transcript that carries a `dbID`
/// and is not stored across the origin of a circular sequence.
///
/// Ensembl VEP's dump writes into a megabase's file every transcript of each
/// gene that overlaps the megabase
/// (`AnnotationSource/Database/Transcript.pm:237, 250`), and for a sequence
/// outside the reference assembly (a patch, an alternate locus) it writes the
/// files of the megabases that sequence's own region lies in and no other
/// (`Pipeline/DumpVEP/CreateDumpJobs.pm:169`,
/// `Pipeline/DumpVEP/Dumper.pm:135-136`). A transcript of a gene that leaves
/// that region therefore reaches megabases the cache holds no file for, where
/// Ensembl VEP reads nothing (`AnnotationSource/Cache/BaseSerialized.pm:255-263`).
///
/// A transcript without a `dbID` is none of that dump's: a cache built from a
/// gene-model file holds a transcript once, in the file of the megabase it
/// starts in, and it is found by its coordinates.
fn filed_short_of_its_span(transcript: &Transcript, spans: &[(u64, u64)]) -> bool {
    if transcript.db_id.is_none() || transcript.circular_length.is_some() {
        return false;
    }
    let mut spans = spans.to_vec();
    spans.sort_unstable();
    // The lowest base of the transcript no file so far covers.
    let mut uncovered = transcript.start;
    for (first, last) in spans {
        if first > uncovered {
            break;
        }
        uncovered = uncovered.max(last.saturating_add(1));
    }
    uncovered <= transcript.end
}

/// Gives each RefSeq model of a merged cache's shard that has a gene symbol and
/// no HGNC id the id of a same-symbol model of the shard. Ensembl VEP's
/// `merge_features` lends the id that way (`AnnotationType/Transcript.pm:290,
/// 305`: the ids are indexed by symbol, the last one seen per symbol, then
/// written onto every model with that symbol) over the features of each
/// annotation buffer, so its `HGNC_ID` on a RefSeq row changes with the buffer's
/// other records; the dumper's own run of the routine lends nothing across the
/// two sources of a merged cache, so its RefSeq models arrive without an id.
/// Lending once per shard, before the de-duplication keeps a boundary
/// transcript's lowest shard, is the rule the native builder follows
/// (`docs/divergences/hgnc-id-copy.md`), and a shard whose models come from one
/// source holds nothing to copy and is left as loaded.
fn copy_hgnc_ids_to_refseq_models(transcripts: &mut [Transcript]) {
    let holds = |source: &SourceCache| {
        transcripts
            .iter()
            .any(|t| t.source_cache.as_ref() == Some(source))
    };
    if !(holds(&SourceCache::Ensembl) && holds(&SourceCache::RefSeq)) {
        return;
    }
    let mut by_symbol: HashMap<Arc<str>, Arc<str>> = HashMap::new();
    for t in transcripts.iter() {
        if let (Some(symbol), Some(hgnc_id)) = (&t.gene_symbol, &t.hgnc_id) {
            by_symbol.insert(symbol.clone(), hgnc_id.clone());
        }
    }
    for t in transcripts
        .iter_mut()
        .filter(|t| t.source_cache == Some(SourceCache::RefSeq) && t.hgnc_id.is_none())
    {
        t.hgnc_id = t
            .gene_symbol
            .as_ref()
            .and_then(|symbol| by_symbol.get(symbol))
            .cloned();
    }
}

/// Refuses a coding transcript whose mapper places the coding sequence where
/// the transcript's own span cannot hold it. A cache without UTR strings has
/// each one padded with `N` to the mapper's `cdna_coding_start - 1` and to the
/// cDNA length less `cdna_coding_end`, so each padding is bounded here by the
/// transcript's genomic span plus the bases its RNA edits insert, the most
/// bases its cDNA can hold, before any record is annotated against it.
///
/// A coding end before the coding start is not refused. Ensembl core takes the
/// coding end of a model with sequence edits from its `_transl_end` attribute
/// as it stands (`Transcript.pm:1042-1044`), so a cache can hold a coding end
/// below 1, and Ensembl VEP reads the model: every cDNA position lies past
/// that end, so none maps to a coding position (`TranscriptMapper.pm:441`).
/// The transcript loads without a coding model (`has_coding_model` of its
/// facts).
fn check_mapper_bounds(transcript: &Transcript) -> Result<(), String> {
    let Some(vefc) = transcript.vefc.as_ref() else {
        return Ok(());
    };
    // No padding is built for a transcript without a coding sequence.
    if vefc.translateable_seq.as_deref().is_none_or(str::is_empty) {
        return Ok(());
    }
    let Some(mapper) = vefc.mapper.as_ref() else {
        return Ok(());
    };
    // An `_rna_edit` attribute (`<start> <end> <bases>` in cDNA coordinates)
    // that inserts bases lengthens the cDNA the mapper measures past the
    // exons' span; an insertion's end lies before its start and replaces
    // nothing.
    let inserted_bases: u64 = rna_edits(transcript)
        .map(|(start, end, bases)| {
            let replaced = if end >= start {
                (end - start).saturating_add(1)
            } else {
                0
            };
            (bases.len() as u64).saturating_sub(replaced)
        })
        .sum();
    let span = transcript
        .end
        .saturating_sub(transcript.start)
        .saturating_add(1)
        .saturating_add(inserted_bases);
    let cdna_len = mapper
        .exon_coord_mapper
        .pairs
        .iter()
        .map(|p| p.from_end)
        .max()
        .unwrap_or(0);
    let five_prime = mapper.cdna_coding_start.saturating_sub(1);
    let three_prime = cdna_len.saturating_sub(mapper.cdna_coding_end);
    if five_prime > span || three_prime > span {
        return Err(format!(
            "transcript {}: the mapper's coding start {} and end {} on a cDNA of {cdna_len} \
             bases leave a UTR longer than the transcript's {span}-base span",
            transcript.stable_id, mapper.cdna_coding_start, mapper.cdna_coding_end
        ));
    }
    Ok(())
}

/// Refuses a transcript whose record names a codon table that is not known
/// (`vep_core::codon::known_table`): every codon under such a table would read
/// `X`, so the consequence of every coding variant on the transcript would be
/// wrong. Ensembl VEP's own translation under one reads each codon as an
/// apostrophe or as nothing (ids 7, 8 and 17 to 20) or, after a warning, under
/// its last table (25 and above, the known table 26 among them). A record
/// without the field is table 1.
fn check_codon_table(transcript: &Transcript) -> Result<(), String> {
    let Some(vefc) = transcript.vefc.as_ref() else {
        return Ok(());
    };
    if vep_core::codon::known_table(vefc.codon_table) {
        return Ok(());
    }
    Err(format!(
        "codon table {} of transcript {} is not a genetic code vep-rs translates with",
        vefc.codon_table, transcript.stable_id
    ))
}

/// The cached peptide of a transcript on genetic code `table`, read under that
/// code. Ensembl's dumper writes it with `Transcript::translate`, whose
/// `Bio::Tools::CodonTable` reads a code it lacks under another table
/// (`vep_core::codon::module_table`). For such a code each residue that other
/// table gives its codon of `coding` is translated again under the
/// transcript's own code, `M` first where the first codon is a start codon of
/// the code; any other residue, a sequence edit or one a dumper already read
/// under the code, is kept. Under every other table the peptide is returned
/// as it is.
fn peptide_under_own_code(peptide: String, coding: &str, table: u8) -> String {
    use vep_core::codon::{is_start_codon, known_table, module_table, translate_sequence};

    let dumped_under = module_table(table);
    if dumped_under == table || !known_table(table) {
        return peptide;
    }
    let coding = coding.as_bytes();
    let residue = |codon: &[u8], first: bool, table: u8| {
        if first && is_start_codon(codon, table) {
            b'M'
        } else {
            // `translate_sequence` reads an ambiguous codon as
            // `Bio::Tools::CodonTable::translate` does, where
            // `translate_codon_with_table` reads `X`.
            translate_sequence(codon, table)
                .first()
                .copied()
                .unwrap_or(b'X')
        }
    };
    let reread: Vec<u8> = peptide
        .bytes()
        .enumerate()
        .map(|(at, cached)| match coding.get(3 * at..3 * at + 3) {
            Some(codon) if cached == residue(codon, at == 0, dumped_under) => {
                residue(codon, at == 0, table)
            }
            _ => cached,
        })
        .collect();
    String::from_utf8(reread).unwrap_or(peptide)
}

/// The `_rna_edit` attributes of `transcript` as `(start, end, bases)`: the
/// span of the unedited cDNA an edit replaces, `end` one before `start` for an
/// insertion, and the bases it puts there, none for a deletion. A value without
/// two leading integers, which no cache writes, is skipped.
fn rna_edits(transcript: &Transcript) -> impl Iterator<Item = (u64, u64, &str)> {
    transcript
        .attributes
        .iter()
        .filter(|attribute| &*attribute.code == "_rna_edit")
        .filter_map(|attribute| {
            let mut fields = attribute.value.split_whitespace();
            let start: u64 = fields.next()?.parse().ok()?;
            let end: u64 = fields.next()?.parse().ok()?;
            Some((start, end, fields.next().unwrap_or("")))
        })
}

/// A stop codon an `_rna_edit` completes, as the cache alone shows it
/// ([`completed_stop_codon`]).
struct CompletedStopCodon {
    /// The codon's bases, upper-cased.
    codon: [u8; 3],
    /// First genomic base the edit replaces, on the forward strand.
    genomic_start: u64,
    /// Last genomic base the edit replaces, on the forward strand.
    genomic_end: u64,
}

/// The stop codon of a transcript whose last exon ends one or two bases into
/// it, or `None` for any other transcript. The cache holds such a codon as an
/// `_rna_edit` that replaces the exon's last one or two bases by the whole
/// codon, the last of a coding sequence that ends with the cDNA, and as a last
/// mapper pair that places the replaced genomic bases on the codon's last
/// bases: Ensembl's `_load_mapper` reads every edit that lengthens the cDNA as
/// an insertion at the edit's first base and compares no base
/// (TranscriptMapper.pm:179, :209-211 and :169 at release/116).
fn completed_stop_codon(transcript: &Transcript) -> Option<CompletedStopCodon> {
    transcript.translation.as_ref()?;
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let [.., before, last] = mapper.exon_coord_mapper.pairs.as_slice() else {
        return None;
    };
    let replaced = last.from_end.checked_sub(last.from_start)?.checked_add(1)?;
    if replaced > 2 || last.to_end.checked_sub(last.to_start)? != replaced - 1 {
        return None;
    }
    // The pair ends the cDNA and the coding sequence, in a whole codon whose
    // first bases no pair holds.
    let coding_start = mapper.cdna_coding_start;
    let codon_start = mapper.cdna_coding_end.checked_sub(2)?;
    let phase = u64::try_from(mapper.start_phase).unwrap_or(0);
    if coding_start == 0
        || codon_start < coding_start
        || ((codon_start - coding_start) % 3 + phase % 3) % 3 != 0
        || last.from_end != mapper.cdna_coding_end
        || before.from_end.checked_add(1) != Some(codon_start)
    {
        return None;
    }
    let ends_transcript = match transcript.strand {
        Strand::Forward => last.ori == 1 && last.to_end == transcript.end,
        Strand::Reverse => last.ori == -1 && last.to_start == transcript.start,
    };
    if !ends_transcript {
        return None;
    }
    let coding = vefc.translateable_seq.as_deref()?.as_bytes();
    let mut codon: [u8; 3] = coding
        .get(coding.len().checked_sub(3)?..)?
        .try_into()
        .ok()?;
    codon.make_ascii_uppercase();
    if !vep_core::codon::is_ter_codon(&codon, vefc.codon_table) {
        return None;
    }
    // An edit's span counts the unedited cDNA: the edits before it move it by
    // the bases they add or remove (`$edit_shift`, TranscriptMapper.pm:179).
    let edits: Vec<(u64, u64, &str)> = rna_edits(transcript).collect();
    let length_change = |&(start, end, bases): &(u64, u64, &str)| {
        bases.len() as i128 - (end as i128 - start as i128 + 1)
    };
    let completes = |&(start, end, bases): &(u64, u64, &str)| {
        let shift: i128 = edits
            .iter()
            .filter(|edit| edit.0 < start)
            .map(length_change)
            .sum();
        end.checked_sub(start) == Some(replaced - 1)
            && bases.as_bytes().eq_ignore_ascii_case(&codon)
            && start as i128 + shift == codon_start as i128
    };
    edits.iter().any(completes).then_some(CompletedStopCodon {
        codon,
        genomic_start: last.to_start,
        genomic_end: last.to_end,
    })
}

/// Pairs the genome's last coding bases with the bases of the stop codon they
/// are, on a transcript whose stop codon an `_rna_edit` completes
/// ([`completed_stop_codon`]), and returns `true` when it moves the pair.
///
/// `genome` reads a span of a chromosome from the run's reference sequence, on
/// the forward strand. Where the stop codon begins with the replaced genomic
/// bases and no other placement of them inside the codon is an identity, the
/// edit appends the codon's remaining bases after the exon's last genomic
/// base: the last pair moves onto the codon's first bases and joins the pair
/// before it. The codon's remaining bases keep their cDNA positions and their
/// place in the coding sequence and have no genomic partner, so the coding end
/// stays the cDNA's last base and the genomic coding bounds stay as loaded. A
/// transcript whose placement the reference sequence does not show keeps the
/// cache's pairs.
pub(crate) fn pair_completed_stop_codon(
    transcript: &mut Transcript,
    genome: impl FnOnce(&str, u64, u64) -> Option<Vec<u8>>,
) -> bool {
    let Some(stop) = completed_stop_codon(transcript) else {
        return false;
    };
    let Some(mut held) = genome(&transcript.chr, stop.genomic_start, stop.genomic_end) else {
        return false;
    };
    if held.len() as u64 != stop.genomic_end - stop.genomic_start + 1 {
        return false;
    }
    if transcript.strand == Strand::Reverse {
        held = vep_core::codon::reverse_complement(&held);
    }
    held.make_ascii_uppercase();
    let mut identities = stop.codon.windows(held.len()).map(|bases| bases == held);
    if identities.next() != Some(true) || identities.any(|identity| identity) {
        return false;
    }
    let Some(mapper) = transcript.vefc.as_mut().and_then(|v| v.mapper.as_mut()) else {
        return false;
    };
    let mut pairs = mapper.exon_coord_mapper.pairs.clone();
    let Some(mut moved) = pairs.pop() else {
        return false;
    };
    let appended = 3 - held.len() as u64;
    moved.from_start -= appended;
    moved.from_end -= appended;
    match pairs.last_mut() {
        Some(before)
            if before.ori == 1
                && moved.ori == 1
                && before.to_end.checked_add(1) == Some(moved.to_start) =>
        {
            before.from_end = moved.from_end;
            before.to_end = moved.to_end;
        }
        Some(before)
            if before.ori == -1
                && moved.ori == -1
                && moved.to_end.checked_add(1) == Some(before.to_start) =>
        {
            before.from_end = moved.from_end;
            before.to_start = moved.to_start;
        }
        _ => pairs.push(moved),
    }
    mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
    true
}

/// [`pair_completed_stop_codon`] on each of `transcripts`, the genome read
/// from `reference`.
pub(crate) fn pair_completed_stop_codons(
    transcripts: &mut [Transcript],
    reference: &vep_fasta::IndexedFasta,
) {
    for transcript in transcripts {
        pair_completed_stop_codon(transcript, |chr, start, end| {
            reference.sequence(chr, start, end)
        });
    }
}

/// A copy of `transcript` with its completed stop codon paired from
/// `reference` ([`pair_completed_stop_codon`]), for a transcript an index
/// already shares; `None` when the pair does not move, as on a transcript a
/// load paired.
pub(crate) fn with_completed_stop_codon_paired(
    transcript: &Transcript,
    reference: &vep_fasta::IndexedFasta,
) -> Option<Transcript> {
    completed_stop_codon(transcript)?;
    let mut paired = transcript.clone();
    pair_completed_stop_codon(&mut paired, |chr, start, end| {
        reference.sequence(chr, start, end)
    })
    .then_some(paired)
}

/// Sets every transcript's `id_rank` to its [`perl_transcript_key`]'s position
/// in the byte order of the slice's keys, counted from 1; transcripts sharing a
/// key share the rank. The key's `_<dbID>` suffix decides the order where one
/// stable id prefixes another: `YIL171W-A_mRNA_34173` precedes
/// `YIL171W_29222` although `YIL171W` alone would precede `YIL171W-A_mRNA`,
/// and the two rank the same way in a cache whose records hold no database
/// id. Rank 0 stays the mark of a transcript no loader ranked, which the
/// writers order by its key instead.
pub(crate) fn assign_id_ranks(transcripts: &mut [Transcript]) {
    let keys: Vec<String> = transcripts.iter().map(perl_transcript_key).collect();
    let mut by_key: Vec<usize> = (0..transcripts.len()).collect();
    by_key.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut rank = 1u32;
    let mut previous: Option<&str> = None;
    for i in by_key {
        if previous.is_some_and(|p| p != keys[i]) {
            rank += 1;
        }
        previous = Some(&keys[i]);
        transcripts[i].id_rank = rank;
    }
}

/// The key Ensembl VEP sorts a variant's transcript variations by
/// (`_get_transcript_key`, BaseVariationFeature.pm:310-321): the stable id,
/// and `_<dbID>` appended when the transcript has a non-zero database id. A
/// transcript of a gene-model file ends in the file's place instead
/// ([`vep_core::transcript::file_order_key_suffix`]). A record without a
/// database id, which a cache built from a GTF or GFF3 holds, ends in `_`
/// alone: Ensembl VEP numbers the transcripts of such a file as it builds
/// them (BaseGXF.pm:489), so every key but the first-built transcript's
/// carries `_` and a number there.
fn perl_transcript_key(transcript: &Transcript) -> String {
    if let Some(order) = transcript
        .source_cache
        .as_ref()
        .and_then(SourceCache::file_order)
    {
        return format!(
            "{}{}",
            transcript.stable_id,
            vep_core::transcript::file_order_key_suffix(order)
        );
    }
    match transcript.db_id {
        Some(0) => transcript.stable_id.to_string(),
        Some(db_id) => format!("{}_{db_id}", transcript.stable_id),
        None => format!("{}_", transcript.stable_id),
    }
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
///
/// The cache's checksum manifest is read from the directory, every file it
/// lists must be present, and every shard is checked against it.
pub fn load_all_transcripts(
    json_cache_dir: &str,
) -> anyhow::Result<HashMap<String, Vec<Transcript>>> {
    load_all_transcripts_in_set(json_cache_dir, None)
}

/// [`load_all_transcripts`] keeping the models of `set` alone, read from each
/// parsed record before its transcript is built ([`load_transcripts_in_set`]).
/// This load covers the whole cache, so a GENCODE set no record of it carries
/// is warned about here ([`GencodeTally::warning`]).
pub(crate) fn load_all_transcripts_in_set(
    json_cache_dir: &str,
    set: Option<&TranscriptSet>,
) -> anyhow::Result<HashMap<String, Vec<Transcript>>> {
    load_all_transcripts_in_set_with(json_cache_dir, set, false)
}

/// [`load_all_transcripts_in_set`] keeping each exon's stable id under
/// `keep_exon_ids` ([`load_transcripts_in_set_with`]).
pub(crate) fn load_all_transcripts_in_set_with(
    json_cache_dir: &str,
    set: Option<&TranscriptSet>,
    keep_exon_ids: bool,
) -> anyhow::Result<HashMap<String, Vec<Transcript>>> {
    let by_chr = load_all_transcripts_without_status_lines(json_cache_dir, set, keep_exon_ids)?;
    write_loaded_lines(&by_chr);
    Ok(by_chr)
}

/// [`load_all_transcripts_in_set_with`] without its status lines, for a caller
/// that adds transcripts to the load before it writes them
/// (`write_loaded_lines`).
pub(crate) fn load_all_transcripts_without_status_lines(
    json_cache_dir: &str,
    set: Option<&TranscriptSet>,
    keep_exon_ids: bool,
) -> anyhow::Result<HashMap<String, Vec<Transcript>>> {
    let shards = enumerate_transcript_shards(json_cache_dir)?;
    let manifest = CacheManifest::load(json_cache_dir).map_err(manifest_error)?;
    manifest.require_complete().map_err(manifest_error)?;
    // Parallel across chromosomes, not only across one chromosome's shards: a
    // serial outer loop drains the pool at every chromosome boundary and runs
    // each tail shard alone. Rayon nests safely here because this eager path
    // runs once at startup outside any `pool.install` scope; the lazy path
    // pre-warms serially for that reason (`LazyTranscriptIndexes::prewarm`).
    let loaded: HashMap<String, LoadedTranscripts> = shards
        .par_iter()
        .map(
            |(chr, paths)| -> anyhow::Result<(String, LoadedTranscripts)> {
                Ok((
                    chr.clone(),
                    load_transcripts_in_set_with(chr, paths, &manifest, set, keep_exon_ids, true)?,
                ))
            },
        )
        .collect::<anyhow::Result<HashMap<_, _>>>()?;
    let mut gencode = GencodeTally::default();
    let mut by_chr: HashMap<String, Vec<Transcript>> = HashMap::with_capacity(loaded.len());
    for (chr, load) in loaded {
        gencode.add(load.gencode);
        by_chr.insert(chr, load.transcripts);
    }
    if let Some(warning) = set
        .and_then(TranscriptSet::gencode_attribute)
        .and_then(|attribute| gencode.warning(attribute))
    {
        tracing::warn!("{warning}");
    }
    Ok(by_chr)
}

/// Writes the status lines of a whole-cache load (`loaded_lines`), which count
/// the transcripts `by_chr` holds when it is called.
pub(crate) fn write_loaded_lines(by_chr: &HashMap<String, Vec<Transcript>>) {
    for line in loaded_lines(by_chr) {
        info!("{line}");
    }
}

/// One chromosome's status line: its name as the cache spells it, with nothing
/// put before it, and the number of transcripts loaded for it.
pub(crate) fn chromosome_line(chromosome: &str, transcripts: usize) -> String {
    format!(
        "  chromosome {chromosome}: {transcripts} transcript{}",
        if transcripts == 1 { "" } else { "s" }
    )
}

/// The status lines of a whole-cache load: the number of chromosomes read, then
/// `chromosome_line` for each in name order.
fn loaded_lines(by_chr: &HashMap<String, Vec<Transcript>>) -> Vec<String> {
    let mut chrs: Vec<&String> = by_chr.keys().collect();
    chrs.sort();
    let mut lines = vec![format!(
        "Loaded transcripts for {} chromosome{}",
        by_chr.len(),
        if by_chr.len() == 1 { "" } else { "s" }
    )];
    lines.extend(
        chrs.into_iter()
            .map(|chr| chromosome_line(chr, by_chr[chr].len())),
    );
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript_index::LazyTranscriptIndexes;
    use crate::variation_source::VariationSource;
    use std::io::Write;
    use vep_core::variation::CachedVariation;
    use vep_io::output::fields::{ClinSigAllele, FieldOptions};

    /// A manifest that checks nothing, for the shards these tests write.
    fn unchecked() -> CacheManifest {
        CacheManifest::absent("")
    }

    /// A directory without `transcripts/`, or whose `transcripts/` holds no
    /// shard, loads nothing and says so: an empty transcript set would annotate
    /// every variant as intergenic and exit 0.
    #[test]
    fn test_load_cache_missing_dir() {
        let error = load_all_transcripts("/nonexistent/path").unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("no cache found at /nonexistent/path: no transcripts directory"),
            "{error}"
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().to_str().unwrap();
        assert!(load_all_transcripts(path).is_err());
        std::fs::create_dir(dir.path().join("transcripts")).unwrap();
        let error = load_all_transcripts(path).unwrap_err();
        assert!(
            error.to_string().contains("no transcript shards under"),
            "{error}"
        );
    }

    /// A copy of a manifest-carrying cache that lacks a listed shard is refused
    /// before any shard is read, naming the file, as a run's cache resolution
    /// refuses it.
    #[test]
    fn eager_load_refuses_a_cache_missing_a_listed_shard() {
        let dir = cache_with_shards(
            "vep_jc_incomplete_",
            &[
                ("21", "[]"),
                ("22", &format!("[{}]", coding_transcript_json(""))),
            ],
        );
        CacheManifest::write(dir.path()).unwrap();
        let path = dir.path().to_str().unwrap();
        assert_eq!(load_all_transcripts(path).unwrap().len(), 2);
        let shard = dir
            .path()
            .join("transcripts")
            .join("22")
            .join("1-1000000.json");
        std::fs::remove_file(&shard).unwrap();
        let error = load_all_transcripts(path).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "ERROR: [cache] transcripts/22/1-1000000.json is listed in CHECKSUMS.sha256 but missing from {path}"
            )
        );
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

    /// The status lines of a whole-cache load name each chromosome as the cache
    /// spells it, in name order, and count one chromosome and one transcript in
    /// the singular.
    #[test]
    fn load_lines_name_a_chromosome_as_the_cache_spells_it() {
        let dir = tempfile::tempdir().unwrap();
        let add = |chr: &str, transcripts: usize| {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            let entries: Vec<String> = (0..transcripts)
                .map(|i| {
                    format!(
                        r#"{{"stable_id":"T{chr}_{i}","gene_stable_id":"G{chr}","chr":"{chr}","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#,
                        start = 100 + i * 10,
                        end = 200 + i * 10
                    )
                })
                .collect();
            std::fs::write(
                chr_dir.join("1-1000000.json"),
                format!("[{}]", entries.join(",")),
            )
            .unwrap();
        };
        let lines = || loaded_lines(&load_all_transcripts(dir.path().to_str().unwrap()).unwrap());
        add("rDNA", 2);
        assert_eq!(
            lines(),
            [
                "Loaded transcripts for 1 chromosome",
                "  chromosome rDNA: 2 transcripts",
            ]
        );
        add("211000022278760", 1);
        add("mitochondrion_genome", 38);
        assert_eq!(
            lines(),
            [
                "Loaded transcripts for 3 chromosomes",
                "  chromosome 211000022278760: 1 transcript",
                "  chromosome mitochondrion_genome: 38 transcripts",
                "  chromosome rDNA: 2 transcripts",
            ]
        );
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

        let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
        let mut ids: Vec<&str> = txs.iter().map(|t| &*t.stable_id).collect();
        ids.sort();
        assert_eq!(ids, vec!["ENST_A", "ENST_NODB", "ENST_STRADDLE"]);
    }

    /// The `source_cache` key a merged cache's region files carry
    /// (`Pipeline/DumpVEP/MergeVEP.pm:103, 110`) loads as the cache it names;
    /// a transcript without the key, every single-source cache's, carries
    /// none; a value that is neither is refused under the `cache` component,
    /// naming the transcript.
    #[test]
    fn source_cache_key_loads_as_ensembl_or_refseq() {
        let refseq = transcript_from(&minimal(r#","source_cache":"RefSeq""#));
        assert_eq!(refseq.source_cache, Some(SourceCache::RefSeq));
        let ensembl = transcript_from(&minimal(r#","source_cache":"Ensembl""#));
        assert_eq!(ensembl.source_cache, Some(SourceCache::Ensembl));
        assert_eq!(transcript_from(&minimal("")).source_cache, None);

        let unknown = minimal(r#","source_cache":"refseq""#);
        let dir = cache_with_shards("vep_jc_source_", &[("21", &format!("[{unknown}]"))]);
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let err = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap_err();
        assert_eq!(
            CliError::from_run(&err),
            "ERROR: [cache] transcript ENST00000400375: source_cache `refseq` is neither \
             Ensembl nor RefSeq"
        );
    }

    /// The `bam_edit_status` `apply_edits` leaves on a model it tried to
    /// correct (`AnnotationType/Transcript.pm:601, 621`) loads as `ok` or
    /// `failed`, and the edited cDNA it stores on an `ok` model
    /// (`:617`) is the transcript's spliced sequence before any assembly,
    /// on a non-coding transcript too, where no sequence could be assembled;
    /// an unknown status is refused naming the transcript.
    #[test]
    fn bam_edit_status_and_cached_spliced_seq_load() {
        let edited = transcript_from(&minimal(
            r#","bam_edit_status":"ok","exons":[{"start":44538981,"end":44538990}],
            "variation_effect_feature_cache":{"spliced_seq":"ACGTACGTAAAC",
            "mapper":{"start_phase":-1,"pairs":[{"from_start":1,"from_end":10,"to_start":44538981,"to_end":44538990,"ori":-1}]}}"#,
        ));
        assert_eq!(edited.bam_edit_status, Some(BamEditStatus::Ok));
        assert!(!edited.facts().has_coding_model);
        assert_eq!(edited.spliced_seq(), Some("ACGTACGTAAAC"));
        assert_eq!(
            edited.vefc.as_ref().unwrap().spliced_seq.as_deref(),
            Some("ACGTACGTAAAC")
        );

        let coding_edited = transcript_from(&minimal(
            r#","bam_edit_status":"ok","translation":{"stable_id":"ENSP1","start":4,"end":12},"exons":[{"start":44538981,"end":44538995}],
            "variation_effect_feature_cache":{"spliced_seq":"GGAATGCCCTAATTTA","five_prime_utr":"GGA","translateable_seq":"ATGCCCTAA","three_prime_utr":"TTT",
            "mapper":{"start_phase":-1,"cdna_coding_start":4,"cdna_coding_end":12,"pairs":[{"from_start":1,"from_end":15,"to_start":44538981,"to_end":44538995,"ori":-1}]}}"#,
        ));
        assert_eq!(coding_edited.spliced_seq(), Some("GGAATGCCCTAATTTA"));

        let failed = transcript_from(&minimal(r#","bam_edit_status":"failed""#));
        assert_eq!(failed.bam_edit_status, Some(BamEditStatus::Failed));
        assert_eq!(failed.spliced_seq(), None);
        assert_eq!(transcript_from(&minimal("")).bam_edit_status, None);

        let jt: JsonTranscript =
            serde_json::from_str(&minimal(r#","bam_edit_status":"OK""#)).unwrap();
        assert_eq!(
            convert_transcript(jt, "21", &mut AttributeCodes::default()).unwrap_err(),
            "transcript ENST00000400375: bam_edit_status `OK` is neither ok nor failed"
        );
    }

    /// A merged cache holds the Ensembl core's dbIDs and the RefSeq
    /// otherfeatures' side by side, so two transcripts sharing a dbID from
    /// different sources are both kept, while the same transcript arriving
    /// from two shards still loads once.
    #[test]
    fn merged_cache_duplicates_are_keyed_by_source_and_db_id() {
        let tx = |id: &str, source: &str, start: u64| {
            format!(
                r#"{{"stable_id":"{id}","dbID":7,"source_cache":"{source}","gene_stable_id":"ENSG","chr":"21","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#,
                end = start + 100
            )
        };
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_merged_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!(
                "[{},{}]",
                tx("ENST_SEVEN", "Ensembl", 100),
                tx("NM_SEVEN", "RefSeq", 999_950)
            ),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("1000001-2000000.json"),
            format!("[{}]", tx("NM_SEVEN", "RefSeq", 999_950)),
        )
        .unwrap();
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
        let mut loaded: Vec<(&str, Option<SourceCache>)> = txs
            .iter()
            .map(|t| (&*t.stable_id, t.source_cache.clone()))
            .collect();
        loaded.sort_by(|a, b| a.0.cmp(b.0));
        assert_eq!(
            loaded,
            vec![
                ("ENST_SEVEN", Some(SourceCache::Ensembl)),
                ("NM_SEVEN", Some(SourceCache::RefSeq)),
            ]
        );
    }

    /// A transcript whose stable id another transcript of the chromosome
    /// shares is held to the spans of the shards it was read from, one span a
    /// shard; a transcript whose id is its own is held to none, also when two
    /// shards hold it, and so is one read from a shard whose name states no
    /// span.
    #[test]
    fn a_shared_stable_id_holds_a_transcript_to_the_shards_it_was_read_from() {
        let tx = |id: &str, db_id: u32, source: &str| {
            format!(
                r#"{{"stable_id":"{id}","dbID":{db_id},"gene_stable_id":"G{db_id}","chr":"21","start":999000,"end":1001000,"strand":1,"biotype":"lncRNA","source":"{source}"}}"#
            )
        };
        let load = |files: &[(&str, String)]| {
            let dir = tempfile::Builder::new()
                .prefix("vep_jc_shared_id_")
                .tempdir()
                .unwrap();
            let chr_dir = dir.path().join("transcripts").join("21");
            std::fs::create_dir_all(&chr_dir).unwrap();
            for (name, records) in files {
                std::fs::write(chr_dir.join(name), format!("[{records}]")).unwrap();
            }
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            let mut held: Vec<_> = load_transcripts_for_chr("21", &shards["21"], &unchecked())
                .unwrap()
                .iter()
                .map(|t| {
                    (
                        t.db_id.unwrap(),
                        t.region_files.as_deref().map(<[_]>::to_vec),
                    )
                })
                .collect();
            held.sort();
            held
        };

        let alignment = tx("XR_SHARED", 11, "ensembl");
        let of_the_gene = tx("XR_SHARED", 22, "Gnomon");
        let alone = tx("XR_ALONE", 33, "Gnomon");
        assert_eq!(
            load(&[
                ("1-1000000.json", format!("{alignment},{alone}")),
                (
                    "1000001-2000000.json",
                    format!("{alignment},{of_the_gene},{alone}")
                ),
            ]),
            [
                (11, Some(vec![(1, 1_000_000), (1_000_001, 2_000_000)])),
                (22, Some(vec![(1_000_001, 2_000_000)])),
                (33, None),
            ]
        );
        assert_eq!(
            load(&[
                ("1-1000000.json", alignment.clone()),
                ("upper.json", of_the_gene.clone()),
            ]),
            [(11, Some(vec![(1, 1_000_000)])), (22, None)]
        );
    }

    /// A transcript with a `dbID` whose shards leave a base of its span
    /// uncovered is held to the spans of those shards, whether it lies partly
    /// or wholly past them; one whose shards cover its span is held to none,
    /// and so are a transcript without a `dbID`, one stored across the origin
    /// of a circular sequence and one read from a shard whose name states no
    /// span.
    #[test]
    fn a_transcript_its_shards_hold_short_of_its_span_is_held_to_them() {
        let tx = |id: &str, db_id: Option<u32>, start: u64, end: u64| {
            let db_id = db_id.map_or(String::new(), |db_id| format!(r#""dbID":{db_id},"#));
            format!(
                r#"{{"stable_id":"{id}",{db_id}"gene_stable_id":"G_{id}","chr":"21","start":{start},"end":{end},"strand":1,"biotype":"lncRNA","source":"ensembl"}}"#
            )
        };
        let load = |files: &[(&str, String)]| {
            let dir = tempfile::Builder::new()
                .prefix("vep_jc_short_of_span_")
                .tempdir()
                .unwrap();
            let chr_dir = dir.path().join("transcripts").join("21");
            std::fs::create_dir_all(&chr_dir).unwrap();
            for (name, records) in files {
                std::fs::write(chr_dir.join(name), format!("[{records}]")).unwrap();
            }
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            let mut held: Vec<_> = load_transcripts_for_chr("21", &shards["21"], &unchecked())
                .unwrap()
                .iter()
                .map(|t| {
                    (
                        t.stable_id.to_string(),
                        t.region_files.as_deref().map(<[_]>::to_vec),
                    )
                })
                .collect();
            held.sort();
            held
        };

        let across_the_end = tx("ACROSS_THE_END", Some(1), 1_999_000, 2_001_000);
        let past_the_end = tx("PAST_THE_END", Some(2), 2_007_000, 2_009_000);
        let across_two_files = tx("ACROSS_TWO_FILES", Some(3), 999_000, 1_001_000);
        let before_the_start = tx("BEFORE_THE_START", Some(4), 999_000, 1_001_000);
        let without_a_database_id = tx("WITHOUT_A_DATABASE_ID", None, 1_999_000, 2_001_000);
        assert_eq!(
            load(&[
                ("1-1000000.json", across_two_files.clone()),
                (
                    "1000001-2000000.json",
                    format!(
                        "{across_the_end},{past_the_end},{across_two_files},\
                         {without_a_database_id}"
                    )
                ),
            ]),
            [
                (
                    "ACROSS_THE_END".to_string(),
                    Some(vec![(1_000_001, 2_000_000)])
                ),
                ("ACROSS_TWO_FILES".to_string(), None),
                (
                    "PAST_THE_END".to_string(),
                    Some(vec![(1_000_001, 2_000_000)])
                ),
                ("WITHOUT_A_DATABASE_ID".to_string(), None),
            ]
        );
        assert_eq!(
            load(&[("1000001-2000000.json", before_the_start.clone())]),
            [(
                "BEFORE_THE_START".to_string(),
                Some(vec![(1_000_001, 2_000_000)])
            )]
        );
        assert_eq!(
            load(&[("upper.json", across_the_end.clone())]),
            [("ACROSS_THE_END".to_string(), None)]
        );
        // The pair across the origin starts at 1999911 and the stored start is
        // -89, so the sequence holds 2,000,000 bases and the transcript loads
        // at 1999911 to 2000060, past the one shard that holds it.
        let across_the_origin = r#"{"stable_id":"ACROSS_THE_ORIGIN","dbID":5,"gene_stable_id":"G_ACROSS_THE_ORIGIN","chr":"21","start":-89,"end":60,"strand":1,"biotype":"lncRNA","source":"ensembl","variation_effect_feature_cache":{"mapper":{"pairs":[{"from_start":1,"from_end":150,"to_start":1999911,"to_end":60,"ori":1}]}}}"#;
        assert_eq!(
            load(&[("1-1000000.json", across_the_origin.to_string())]),
            [("ACROSS_THE_ORIGIN".to_string(), None)]
        );
    }

    /// A merged cache's RefSeq model with a gene symbol and no HGNC id takes
    /// the id of a same-symbol model of its own shard; a RefSeq model whose
    /// symbol no model of the shard carries with an id, or without a symbol,
    /// stays without one, and an Ensembl model is left as loaded. A cache whose
    /// models come from one source copies nothing. A transcript stored in two
    /// shards prints its lowest shard's copy, so a donor in the lower shard
    /// alone lends the id and a donor in the upper shard alone does not.
    #[test]
    fn merged_cache_refseq_models_take_the_hgnc_id_of_a_same_symbol_model() {
        let tx = |id: &str, db_id: u32, source: &str, gene: &str| {
            format!(
                r#"{{"stable_id":"{id}","dbID":{db_id},"source_cache":"{source}","gene_stable_id":"G{db_id}",{gene}"chr":"21","start":100,"end":200,"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#
            )
        };
        let named =
            r#""gene_symbol":"KRTAP10-1","gene_symbol_source":"HGNC","gene_hgnc_id":"HGNC:22966","#;
        let symbol_only = r#""gene_symbol":"KRTAP10-1","#;
        let other_symbol = r#""gene_symbol":"LOC102724159","#;
        let other_named =
            r#""gene_symbol":"MRPS6","gene_symbol_source":"HGNC","gene_hgnc_id":"HGNC:14051","#;
        let other_symbol_only = r#""gene_symbol":"MRPS6","#;
        let load = |prefix: &str, shard_records: &[&[String]]| {
            let dir = tempfile::Builder::new().prefix(prefix).tempdir().unwrap();
            let chr_dir = dir.path().join("transcripts").join("21");
            std::fs::create_dir_all(&chr_dir).unwrap();
            for (i, records) in shard_records.iter().enumerate() {
                let file_name = format!("{}-{}.json", i * 1_000_000 + 1, (i + 1) * 1_000_000);
                std::fs::write(chr_dir.join(file_name), format!("[{}]", records.join(",")))
                    .unwrap();
            }
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
            let mut ids: Vec<(String, Option<String>)> = txs
                .iter()
                .map(|t| {
                    (
                        t.stable_id.to_string(),
                        t.hgnc_id.as_deref().map(str::to_string),
                    )
                })
                .collect();
            ids.sort();
            ids
        };
        let merged = load(
            "vep_jc_hgnc_merged_",
            &[&[
                tx("ENST_NAMED", 1, "Ensembl", named),
                tx("ENST_UNNAMED", 2, "Ensembl", symbol_only),
                tx("NM_SAME_SYMBOL", 3, "RefSeq", symbol_only),
                tx("NM_OTHER_SYMBOL", 4, "RefSeq", other_symbol),
                tx("NM_NO_SYMBOL", 5, "RefSeq", ""),
            ]],
        );
        assert_eq!(
            merged,
            vec![
                ("ENST_NAMED".to_string(), Some("HGNC:22966".to_string())),
                ("ENST_UNNAMED".to_string(), None),
                ("NM_NO_SYMBOL".to_string(), None),
                ("NM_OTHER_SYMBOL".to_string(), None),
                ("NM_SAME_SYMBOL".to_string(), Some("HGNC:22966".to_string())),
            ]
        );
        let refseq_only = load(
            "vep_jc_hgnc_refseq_",
            &[&[
                tx("NM_NAMED", 1, "RefSeq", named),
                tx("NM_SAME_SYMBOL", 2, "RefSeq", symbol_only),
            ]],
        );
        assert_eq!(
            refseq_only,
            vec![
                ("NM_NAMED".to_string(), Some("HGNC:22966".to_string())),
                ("NM_SAME_SYMBOL".to_string(), None),
            ]
        );
        // The donor lends within its shard alone: NM_LOW and NM_HIGH lie in both
        // shards and take the lower shard's copy, NM_FAR lies in the upper shard
        // with the donor in the lower one, and NM_UPPER lies in both shards with
        // its donor in the upper one alone.
        let two_shards = load(
            "vep_jc_hgnc_shards_",
            &[
                &[
                    tx("ENST_NAMED", 1, "Ensembl", named),
                    tx("NM_LOW", 2, "RefSeq", symbol_only),
                    tx("NM_HIGH", 3, "RefSeq", symbol_only),
                    tx("NM_UPPER", 7, "RefSeq", other_symbol_only),
                ],
                &[
                    tx("ENST_UNNAMED", 4, "Ensembl", symbol_only),
                    tx("NM_LOW", 2, "RefSeq", symbol_only),
                    tx("NM_HIGH", 3, "RefSeq", symbol_only),
                    tx("NM_FAR", 5, "RefSeq", symbol_only),
                    tx("ENST_OTHER", 6, "Ensembl", other_named),
                    tx("NM_UPPER", 7, "RefSeq", other_symbol_only),
                ],
            ],
        );
        assert_eq!(
            two_shards,
            vec![
                ("ENST_NAMED".to_string(), Some("HGNC:22966".to_string())),
                ("ENST_OTHER".to_string(), Some("HGNC:14051".to_string())),
                ("ENST_UNNAMED".to_string(), None),
                ("NM_FAR".to_string(), None),
                ("NM_HIGH".to_string(), Some("HGNC:22966".to_string())),
                ("NM_LOW".to_string(), Some("HGNC:22966".to_string())),
                ("NM_UPPER".to_string(), None),
            ]
        );
    }

    /// The set read from the parsed record keeps the transcripts the set read
    /// from the built transcript keeps, in the same order with the same
    /// relative ranks, on a merged shard: the RefSeq set, the merged set, the
    /// merged set with `--all_refseq` and the GENCODE basic set on a merged
    /// cache. A dropped record is never converted, so a dropped record whose
    /// mapper the loader would refuse loads under the set that drops it and
    /// refuses the cache without one; a `source_cache` value the cache cannot
    /// hold refuses it either way.
    #[test]
    fn set_read_before_conversion_keeps_what_the_built_transcript_set_keeps() {
        use crate::cache_resolve::SourceType;
        use crate::transcript_index::TranscriptSet;

        let tx = |id: &str, source: &str, start: u64, attributes: &str, extra: &str| {
            format!(
                r#"{{"stable_id":"{id}","dbID":{start},"source_cache":"{source}","gene_stable_id":"ENSG","chr":"21","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl","attributes":[{attributes}]{extra}}}"#,
                end = start + 100
            )
        };
        let basic = r#"{"code":"gencode_basic","value":"GENCODE basic"}"#;
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_set_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!(
                "[{},{},{},{},{}]",
                tx("ENST_A", "Ensembl", 100, basic, ""),
                tx("ENST_B", "Ensembl", 200, "", ""),
                tx("NM_003", "RefSeq", 300, basic, ""),
                tx("LOC_D", "RefSeq", 400, basic, ""),
                tx("NR_005", "RefSeq", 500, "", ""),
            ),
        )
        .unwrap();
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let ids = |txs: &[Transcript]| -> Vec<(String, u32)> {
            txs.iter()
                .map(|t| (t.stable_id.to_string(), t.id_rank))
                .collect()
        };
        let sets = [
            TranscriptSet {
                source_type: SourceType::RefSeq,
                ..TranscriptSet::default()
            },
            TranscriptSet {
                source_type: SourceType::Merged,
                ..TranscriptSet::default()
            },
            TranscriptSet {
                source_type: SourceType::Merged,
                all_refseq: true,
                ..TranscriptSet::default()
            },
            TranscriptSet {
                source_type: SourceType::Merged,
                gencode_basic: true,
                ..TranscriptSet::default()
            },
            TranscriptSet::default(),
        ];
        for set in sets {
            let before = load_transcripts_in_set("21", &shards["21"], &unchecked(), Some(&set))
                .unwrap()
                .transcripts;
            let mut after = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
            after.retain(|t| set.accepts(t));
            let before_ids: Vec<String> = ids(&before).into_iter().map(|(id, _)| id).collect();
            let after_ids: Vec<String> = ids(&after).into_iter().map(|(id, _)| id).collect();
            assert_eq!(before_ids, after_ids, "{set:?}");
            // The ranks of the kept transcripts order them alike either way.
            let order = |txs: &[Transcript]| -> Vec<usize> {
                let ranks = ids(txs);
                let mut by_rank: Vec<usize> = (0..ranks.len()).collect();
                by_rank.sort_by_key(|&i| ranks[i].1);
                by_rank
            };
            assert_eq!(order(&before), order(&after), "{set:?}");
        }
        let refseq = load_transcripts_in_set(
            "21",
            &shards["21"],
            &unchecked(),
            Some(&TranscriptSet {
                source_type: SourceType::RefSeq,
                ..TranscriptSet::default()
            }),
        )
        .unwrap()
        .transcripts;
        assert_eq!(
            ids(&refseq),
            vec![("NM_003".to_string(), 1), ("NR_005".to_string(), 2)]
        );

        // A dropped record is never converted: its defects pass unseen under
        // the set that drops it and refuse the cache without one.
        let defective = dir.path().join("transcripts").join("22");
        std::fs::create_dir_all(&defective).unwrap();
        let mapper = r#","variation_effect_feature_cache":{"translateable_seq":"ATG","mapper":{"start_phase":0,"cdna_coding_start":5000,"cdna_coding_end":5002,"pairs":[{"from_start":1,"from_end":6000,"to_start":100,"to_end":200,"ori":1}]}}"#;
        std::fs::write(
            defective.join("1-1000000.json"),
            format!(
                "[{},{}]",
                tx("NM_001", "RefSeq", 100, "", ""),
                tx("LOC_BAD", "RefSeq", 300, "", mapper),
            ),
        )
        .unwrap();
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let refseq_set = TranscriptSet {
            source_type: SourceType::RefSeq,
            ..TranscriptSet::default()
        };
        let kept = load_transcripts_in_set("22", &shards["22"], &unchecked(), Some(&refseq_set))
            .unwrap()
            .transcripts;
        assert_eq!(ids(&kept), vec![("NM_001".to_string(), 1)]);
        let refused = load_transcripts_for_chr("22", &shards["22"], &unchecked())
            .unwrap_err()
            .to_string();
        assert!(refused.contains("transcript LOC_BAD"), "{refused}");
        std::fs::write(
            defective.join("1-1000000.json"),
            format!("[{}]", tx("LOC_BAD", "Other", 300, "", "")),
        )
        .unwrap();
        let refused = load_transcripts_in_set("22", &shards["22"], &unchecked(), Some(&refseq_set))
            .unwrap_err()
            .to_string();
        assert!(
            refused.contains("source_cache `Other` is neither Ensembl nor RefSeq"),
            "{refused}"
        );
    }

    /// Under a GENCODE flag a load counts the records it parsed, over every
    /// shard, and those carrying the set's attribute, a dropped record included,
    /// and owes the warning when none carried it and some was parsed; without
    /// the flag it counts nothing, and no warning is owed for an empty load.
    #[test]
    fn gencode_tally_counts_the_sets_attribute_over_every_shard() {
        use crate::transcript_index::TranscriptSet;

        let tx = |id: &str, start: u64, attributes: &str| {
            format!(
                r#"{{"stable_id":"{id}","gene_stable_id":"ENSG","chr":"21","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl","attributes":[{attributes}]}}"#,
                end = start + 100
            )
        };
        let basic = r#"{"code":"gencode_basic","value":"GENCODE basic"}"#;
        let primary = r#"{"code":"gencode_primary","value":"GENCODE primary"}"#;
        let dir = cache_with_shards(
            "vep_jc_tally_",
            &[
                (
                    "21",
                    &format!("[{},{}]", tx("ENST_A", 100, basic), tx("ENST_B", 200, "")),
                ),
                ("22", &format!("[{}]", tx("ENST_C", 300, primary))),
            ],
        );
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let both: Vec<PathBuf> = shards["21"].iter().chain(&shards["22"]).cloned().collect();
        let basic_set = TranscriptSet {
            gencode_basic: true,
            ..TranscriptSet::default()
        };
        let primary_set = TranscriptSet {
            gencode_primary: true,
            ..TranscriptSet::default()
        };

        let load = |set: Option<&TranscriptSet>| {
            load_transcripts_in_set("21", &both, &unchecked(), set).unwrap()
        };
        let no_flag = load(None);
        assert_eq!(no_flag.transcripts.len(), 3);
        assert_eq!(no_flag.gencode, GencodeTally::default());
        assert_eq!(no_flag.gencode.warning("gencode_basic"), None);

        let basic_load = load(Some(&basic_set));
        assert_eq!(
            basic_load.gencode,
            GencodeTally {
                records: 3,
                carriers: 1
            }
        );
        assert_eq!(basic_load.transcripts.len(), 1);
        assert_eq!(basic_load.gencode.warning("gencode_basic"), None);

        let primary_load = load(Some(&primary_set));
        assert_eq!(
            primary_load.gencode,
            GencodeTally {
                records: 3,
                carriers: 1
            }
        );

        // The chromosome 21 shard alone carries no `gencode_primary`.
        let none =
            load_transcripts_in_set("21", &shards["21"], &unchecked(), Some(&primary_set)).unwrap();
        assert!(none.transcripts.is_empty());
        assert_eq!(
            none.gencode.warning("gencode_primary").as_deref(),
            Some(
                "no transcript of the cache carries the gencode_primary attribute; \
                 every transcript is filtered"
            )
        );
        let mut summed = none.gencode;
        summed.add(primary_load.gencode);
        assert_eq!(
            summed,
            GencodeTally {
                records: 5,
                carriers: 1
            }
        );
        assert_eq!(summed.warning("gencode_primary"), None);
        assert_eq!(GencodeTally::default().warning("gencode_basic"), None);
    }

    /// `id_rank` orders a chromosome's transcripts as Ensembl VEP's transcript
    /// keys order as bytes (the stable id, then `_` and the dbID), whatever order
    /// the shards listed them in, counted from 1 so that no ranked transcript
    /// carries the unranked 0; a transcript written into two shards loads once,
    /// and two transcripts with one id and distinct dbIDs rank by dbID.
    #[test]
    fn loaded_transcripts_rank_by_ensembl_veps_transcript_key() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_rank_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let tx = |id: &str, db_id: u32, start: u64| {
            format!(
                r#"{{"stable_id":"{id}","dbID":"{db_id}","gene_stable_id":"ENSG","chr":"21","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#,
                end = start + 100
            )
        };
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!(
                "[{},{},{},{}]",
                tx("ENST00000300", 1, 100),
                tx("ENST00000010", 2, 200),
                tx("ENST00000200", 4, 300),
                tx("ENST00000200", 3, 999_950)
            ),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("1000001-2000000.json"),
            format!(
                "[{},{}]",
                tx("ENST00000200", 3, 999_950),
                tx("ENST0000002", 5, 1000200)
            ),
        )
        .unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
        let mut ranked: Vec<(&str, Option<u64>, u32)> = txs
            .iter()
            .map(|t| (&*t.stable_id, t.db_id, t.id_rank))
            .collect();
        ranked.sort();
        // Byte order, not numeric: the shorter `ENST0000002` sorts between
        // `ENST00000010` and `ENST00000200`; the `ENST00000200` record with dbID
        // 3 sits in both shards and loads once.
        assert_eq!(
            ranked,
            vec![
                ("ENST00000010", Some(2), 1),
                ("ENST0000002", Some(5), 2),
                ("ENST00000200", Some(3), 3),
                ("ENST00000200", Some(4), 4),
                ("ENST00000300", Some(1), 5),
            ]
        );
        assert!(txs.iter().all(|t| t.id_rank != 0));
        for a in &txs {
            for b in &txs {
                assert_eq!(
                    a.id_rank.cmp(&b.id_rank),
                    perl_transcript_key(a).cmp(&perl_transcript_key(b)),
                    "{} against {}",
                    a.stable_id,
                    b.stable_id
                );
            }
        }
    }

    /// Where one stable id is a prefix of another, the `_<dbID>` suffix decides
    /// the order as it does in Ensembl VEP: `YIL171W-A_mRNA_34173` precedes
    /// `YIL171W_29222` because `-` sorts before `_`, and `YML045W-A_35661`
    /// precedes `YML045W_35617`. A record without a dbID, as a cache built
    /// from a GTF or GFF3 holds it, keys as its id and `_`, so such a pair
    /// ranks the same way; a dbID of zero leaves the key bare
    /// (BaseVariationFeature.pm:316).
    #[test]
    fn a_stable_id_that_prefixes_another_ranks_by_the_dbid_suffix() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_rank_prefix_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("IX");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let tx = |id: &str, db_id: Option<u32>, start: u64| {
            let db_id = match db_id {
                Some(db_id) => format!(r#""dbID":"{db_id}","#),
                None => String::new(),
            };
            format!(
                r#"{{"stable_id":"{id}",{db_id}"gene_stable_id":"G","chr":"IX","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"sgd"}}"#,
                end = start + 100
            )
        };
        std::fs::write(
            chr_dir.join("1-1000000.json"),
            format!(
                "[{},{},{},{},{},{},{},{},{}]",
                tx("YIL171W", Some(29222), 19515),
                tx("YIL171W-A_mRNA", Some(34173), 18260),
                tx("YML045W", Some(35617), 184461),
                tx("YML045W-A", Some(35661), 184461),
                tx("YML044C", None, 180000),
                tx("MODEL1", None, 200000),
                tx("MODEL1-A", None, 200500),
                tx("MODEL2", Some(0), 210000),
                tx("MODEL2-A", Some(8), 210500)
            ),
        )
        .unwrap();

        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        let txs = load_transcripts_for_chr("IX", &shards["IX"], &unchecked()).unwrap();
        let transcript = |id: &str| txs.iter().find(|t| &*t.stable_id == id).unwrap();
        let rank = |id: &str| transcript(id).id_rank;
        assert!(rank("YIL171W-A_mRNA") < rank("YIL171W"));
        assert!(rank("YML045W-A") < rank("YML045W"));
        assert!(rank("YML044C") < rank("YML045W-A"));
        assert_eq!(perl_transcript_key(transcript("YML044C")), "YML044C_");
        assert!(rank("MODEL1-A") < rank("MODEL1"));
        assert_eq!(perl_transcript_key(transcript("MODEL2")), "MODEL2");
        assert!(rank("MODEL2") < rank("MODEL2-A"));
    }

    /// The key the loader ranks a transcript by and the key a row without a
    /// rank is ordered by (`TranscriptConsequence::ensembl_vep_key`) are the
    /// same bytes: with a database id, without one, with a zero one, and for a
    /// transcript of a gene-model file, whatever database id that one carries.
    #[test]
    fn the_rank_key_and_the_row_key_of_a_transcript_are_the_same_bytes() {
        use vep_core::consequence::TranscriptConsequence;
        let file = |order: u32| SourceCache::GeneModelFile {
            short_name: "genes.gtf.gz".into(),
            order,
        };
        let cases: [(Option<u64>, Option<SourceCache>, &str); 7] = [
            (Some(34173), None, "MODEL1_34173"),
            (None, None, "MODEL1_"),
            (Some(0), None, "MODEL1"),
            (Some(7), Some(SourceCache::RefSeq), "MODEL1_7"),
            (None, Some(file(1)), "MODEL1_:0000000001"),
            (Some(7), Some(file(12)), "MODEL1_:0000000012"),
            (Some(0), Some(file(2)), "MODEL1_:0000000002"),
        ];
        for (db_id, source_cache, key) in cases {
            let mut transcript = Transcript::default();
            transcript.stable_id = "MODEL1".into();
            transcript.db_id = db_id;
            transcript.source_cache = source_cache;
            assert_eq!(perl_transcript_key(&transcript), key);
            let row = TranscriptConsequence::for_transcript(&Arc::new(transcript));
            let row_key = String::from_utf8(row.ensembl_vep_key().collect()).unwrap();
            assert_eq!(row_key, key);
        }
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
        let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
        assert_eq!(txs.len(), 1);
        assert_eq!(&*txs[0].stable_id, "ENST_GZ");
    }

    /// A shard that wraps its transcripts in a map of arrays with null
    /// placeholders loads beside an array shard, nulls dropped, keys in sorted
    /// order.
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
        let txs = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap();
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
        let err = load_transcripts_for_chr("21", &shards["21"], &unchecked()).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("1-1000000.json"), "{msg}");
        assert!(msg.contains("EOF"), "{msg}");
        assert!(!msg.contains("expected a sequence"), "{msg}");
    }

    /// A coding end the cache holds as a negative number, as a JSON number or
    /// as a string holding one, reads as 0 on the transcript and in its
    /// mapper, and the coding transcript that carries it loads without a
    /// coding model and without coding bounds on the genome. Ensembl VEP reads
    /// such a model, and no coding position maps to it (ensembl core
    /// `TranscriptMapper.pm:441`: every cDNA position lies past the coding
    /// end).
    #[test]
    fn a_negative_coding_end_loads_without_a_coding_model() {
        for stored in [r#""-20""#, "-20"] {
            let json = coding_transcript_json("")
                .replace(
                    r#""cdna_coding_end":30,"#,
                    &format!(r#""cdna_coding_end":{stored},"#),
                )
                .replace(
                    r#""stable_id":"ENST1","#,
                    &format!(
                        r#""stable_id":"ENST1","cdna_coding_start":1,"cdna_coding_end":{stored},"#
                    ),
                );
            let dir = cache_with_shards("vep_jc_negative_end_", &[("21", &format!("[{json}]"))]);
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            let loaded = load_transcripts_for_chr("21", &shards["21"], &unchecked())
                .unwrap_or_else(|error| panic!("{stored}: {error}"));
            assert_eq!(loaded.len(), 1, "{stored}");
            let transcript = &loaded[0];
            assert_eq!(
                (transcript.cdna_coding_start, transcript.cdna_coding_end),
                (Some(1), Some(0)),
                "{stored}"
            );
            let mapper = transcript.vefc.as_ref().unwrap().mapper.as_ref().unwrap();
            assert_eq!(
                (mapper.cdna_coding_start, mapper.cdna_coding_end),
                (1, 0),
                "{stored}"
            );
            assert!(transcript.translation.is_some(), "{stored}");
            assert!(!transcript.facts().has_coding_model, "{stored}");
            assert_eq!(
                (transcript.coding_region_start, transcript.coding_region_end),
                (None, None),
                "{stored}"
            );
        }
    }

    /// A shard whose mapper places the coding start, or the coding end, a UTR
    /// longer than the transcript's span away is refused under the `cache`
    /// component, naming the transcript, when its chromosome loads; the
    /// 30-base transcript with its coding sequence over the whole span loads,
    /// as do a non-coding transcript whose mapper runs past its span (no
    /// padding is built for it) and a coding one whose cDNA runs past its span
    /// by the bases an RNA edit inserts; an edit whose replaced range fills the
    /// cDNA numbering inserts nothing and is read without a failure.
    #[test]
    fn mapper_utr_longer_than_the_transcript_span_is_refused() {
        let load = |replacements: &[(&str, &str)]| {
            let mut json = coding_transcript_json("");
            for (from, to) in replacements {
                assert!(json.contains(from), "{from}");
                json = json.replace(from, to);
            }
            let dir = cache_with_shards("vep_jc_mapper_", &[("21", &format!("[{json}]"))]);
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            load_transcripts_for_chr("21", &shards["21"], &unchecked()).map(|txs| txs.len())
        };
        // The fixture's one mapper pair ends at cDNA 30, so its `from_end` is
        // the cDNA length.
        let cdna_far_past_span = (r#""from_end":30,"#, r#""from_end":4000000000,"#);
        let cdna_five_bases_past_span = (r#""from_end":30,"#, r#""from_end":65,"#);
        assert_eq!(load(&[]).unwrap(), 1);

        let far_start = load(&[(
            r#""cdna_coding_start":1,"#,
            &format!(r#""cdna_coding_start":{},"#, u64::MAX),
        )])
        .unwrap_err();
        assert_eq!(
            CliError::from_run(&far_start),
            format!(
                "ERROR: [cache] transcript ENST1: the mapper's coding start {} and end 30 on a \
                 cDNA of 30 bases leave a UTR longer than the transcript's 30-base span",
                u64::MAX
            )
        );
        let far_end = load(&[cdna_far_past_span]).unwrap_err();
        assert_eq!(
            CliError::from_run(&far_end),
            "ERROR: [cache] transcript ENST1: the mapper's coding start 1 and end 30 on a cDNA \
             of 4000000000 bases leave a UTR longer than the transcript's 30-base span"
        );

        // A five-base RNA edit insertion lets the cDNA run five bases past the span.
        assert!(load(&[cdna_five_bases_past_span]).is_err());
        let inserting_edit = (
            r#""stable_id":"ENST1","#,
            r#""stable_id":"ENST1","attributes":[{"code":"_rna_edit","value":"5 4 AAAAA"}],"#,
        );
        assert_eq!(
            load(&[cdna_five_bases_past_span, inserting_edit]).unwrap(),
            1
        );
        // An edit replacing cDNA 0 through `u64::MAX` with one base inserts
        // nothing: the transcript loads with its cDNA inside the span and is
        // refused, not abandoned, with the cDNA five bases past it.
        let whole_range = format!(
            r#""stable_id":"ENST1","attributes":[{{"code":"_rna_edit","value":"0 {} A"}}],"#,
            u64::MAX
        );
        let whole_range_edit = (r#""stable_id":"ENST1","#, whole_range.as_str());
        assert_eq!(load(&[whole_range_edit]).unwrap(), 1);
        assert!(load(&[cdna_five_bases_past_span, whole_range_edit]).is_err());

        // The mapper the far-end case refused loads on a transcript without a
        // coding sequence: no padding is built for it, so the bound does not apply.
        let non_coding = (
            r#""translateable_seq":"ATGGCTGCTGCTGCTGCTGCTGCTGCTTGA","#,
            r#""translateable_seq":"","#,
        );
        assert_eq!(load(&[cdna_far_past_span, non_coding]).unwrap(), 1);

        // A coding end before the coding start loads, on a coding transcript
        // as on one without a coding sequence.
        let inverted = (r#""cdna_coding_end":30,"#, r#""cdna_coding_end":0,"#);
        assert_eq!(load(&[inverted]).unwrap(), 1);
        assert_eq!(load(&[inverted, non_coding]).unwrap(), 1);
    }

    /// A shard whose record names a codon table that is not known is refused
    /// under the `cache` component when its chromosome loads, naming the table
    /// and the transcript. BioPerl 1.006924 holds no table at 7, 8 and 17 to
    /// 20 and none from 25 on, and of those 26 alone is known. A record naming
    /// a known table, as a number or a string, loads, as does one without the
    /// field, which is table 1.
    #[test]
    fn a_codon_table_that_is_not_known_is_refused_on_load() {
        let load = |codon_table: Option<&str>| {
            let field = r#""codon_table":1,"#;
            let json = coding_transcript_json("");
            assert!(json.contains(field));
            let json = match codon_table {
                Some(value) => json.replace(field, &format!(r#""codon_table":{value},"#)),
                None => json.replace(field, ""),
            };
            let dir = cache_with_shards("vep_jc_codon_table_", &[("21", &format!("[{json}]"))]);
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            load_transcripts_for_chr("21", &shards["21"], &unchecked())
        };
        for refused in [
            "7", "8", "17", "18", "19", "20", "25", "27", "255", r#""7""#,
        ] {
            let err = load(Some(refused)).unwrap_err();
            assert_eq!(
                CliError::from_run(&err),
                format!(
                    "ERROR: [cache] codon table {} of transcript ENST1 is not a genetic code \
                     vep-rs translates with",
                    refused.trim_matches('"')
                )
            );
        }
        let known = [
            "0", "1", "2", "3", "5", "11", "24", "26", r#""2""#, r#""11""#, r#""26""#,
        ];
        for shipped in known {
            let txs = load(Some(shipped)).unwrap();
            assert_eq!(txs.len(), 1, "{shipped}");
            let want: u8 = shipped.trim_matches('"').parse().unwrap();
            assert_eq!(txs[0].vefc.as_ref().unwrap().codon_table, want, "{shipped}");
        }
        let without_field = load(None).unwrap();
        assert_eq!(without_field[0].vefc.as_ref().unwrap().codon_table, 1);
    }

    /// The cached peptide of a transcript on genetic code 26 is read under
    /// that code. `Bio::Tools::CodonTable` holds no table 26 and reads the id
    /// under its last table, the standard translation with `ATG` its one start
    /// codon (`CodonTable.pm:328-338` and `:412-460`), so a dumper on that
    /// module writes leucine for `CTG` and no `M` over a first `CTG`. Each
    /// residue that table gives its codon is read again under code 26, which
    /// reads `CTG` as alanine and as a start codon; a residue it does not
    /// give, a sequence edit or one already read under code 26, is kept, as is
    /// every residue past the whole codons. A peptide on a table the module
    /// holds, or on an unknown one, is returned as it is.
    #[test]
    fn the_cached_peptide_on_genetic_code_26_is_read_under_that_code() {
        let reread = |peptide: &str, coding: &str, table: u8| {
            peptide_under_own_code(peptide.to_string(), coding, table)
        };
        // ATG GCT CTG CTA GAA TAA
        let coding = "ATGGCTCTGCTAGAATAA";
        assert_eq!(reread("MALLE", coding, 26), "MAALE");
        assert_eq!(reread("MALLE", &coding.to_lowercase(), 26), "MAALE");
        assert_eq!(reread("MAALE", coding, 26), "MAALE");
        assert_eq!(reread("MAULE", coding, 26), "MAULE");
        assert_eq!(reread("MALLEKK", coding, 26), "MAALEKK");
        assert_eq!(reread("", coding, 26), "");
        // CTG GCT TAA
        assert_eq!(reread("LA", "CTGGCTTAA", 26), "MA");
        assert_eq!(reread("MA", "CTGGCTTAA", 26), "MA");
        // ATG CTN YTG TAA: every reading of `CTN` and of `YTG` is leucine under
        // the module's last table, and `CTG` is alanine under code 26.
        assert_eq!(reread("MLL", "ATGCTNYTGTAA", 26), "MXX");
        for table in [1, 4, 11, 12, 24, 25, 27] {
            assert_eq!(reread("MALLE", coding, table), "MALLE", "table {table}");
        }
    }

    /// A cache whose sequence regions name different genetic codes keeps each
    /// transcript on its own code. Ensembl VEP stores the table on every
    /// transcript's record, from its sequence region's `codon_table`
    /// attribute (`AnnotationSource/Database/Transcript.pm:359-363`), and
    /// translates under the record's value
    /// (`BaseTranscriptVariation.pm:1146-1164`), so a region's shard loads
    /// with the table its records name whatever another region's shard names,
    /// before and after that one loads, and a record without the field beside
    /// them is table 1. The pairs are the standard code beside each other
    /// genetic code the ready-made Ensembl VEP caches name, and two
    /// non-standard codes side by side. Each non-standard code but table 11
    /// translates some codon unlike table 1, so a record read under the
    /// default table is a wrong reading; table 11 translates all 64 codons as
    /// table 1 does and marks other start codons.
    #[test]
    fn two_sequences_of_one_cache_keep_their_own_codes() {
        use vep_core::codon::{is_start_codon, translate_codon_with_table};

        // The 64 codons in the order of `Bio::Tools::CodonTable`'s strings.
        let bases = *b"TCAG";
        let mut codons = Vec::with_capacity(64);
        for first in bases {
            for second in bases {
                for third in bases {
                    codons.push([first, second, third]);
                }
            }
        }
        let translates_as_standard = |table: u8| {
            codons.iter().all(|codon| {
                translate_codon_with_table(codon, table) == translate_codon_with_table(codon, 1)
            })
        };
        let starts_as_standard = |table: u8| {
            codons
                .iter()
                .all(|codon| is_start_codon(codon, table) == is_start_codon(codon, 1))
        };

        let record = |region: &str, id: &str, codon_table: Option<u8>| {
            let (stable_id, chr) = (r#""stable_id":"ENST1""#, r#""chr":"21""#);
            let field = r#""codon_table":1,"#;
            let json = coding_transcript_json("");
            assert!(json.contains(stable_id) && json.contains(chr) && json.contains(field));
            let json = json
                .replace(stable_id, &format!(r#""stable_id":"{id}""#))
                .replace(chr, &format!(r#""chr":"{region}""#));
            match codon_table {
                Some(table) => json.replace(field, &format!(r#""codon_table":{table},"#)),
                None => json.replace(field, ""),
            }
        };
        let shard = |region: &str, table: u8| {
            let named = record(region, "named", Some(table));
            let unnamed = record(region, "unnamed", None);
            format!("[{named},{unnamed}]")
        };
        let loaded_as = |region: &str, table: u8| {
            vec![
                (region.to_string(), "named".to_string(), table),
                (region.to_string(), "unnamed".to_string(), 1),
            ]
        };

        assert!(!starts_as_standard(11));
        let beside_standard = [2u8, 3, 4, 5, 6, 9, 11, 12, 13].map(|other| (1u8, other));
        let both_non_standard = [(11, 4), (12, 4)];
        for (first, second) in beside_standard.into_iter().chain(both_non_standard) {
            assert_ne!(first, second);
            for table in [first, second].into_iter().filter(|table| *table != 1) {
                assert_eq!(translates_as_standard(table), table == 11, "table {table}");
            }

            let dir = cache_with_shards(
                "vep_jc_two_codes_",
                &[("1", &shard("1", first)), ("2", &shard("2", second))],
            );
            let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
            let tables = |region: &str| -> Vec<(String, String, u8)> {
                load_transcripts_for_chr(region, &shards[region], &unchecked())
                    .unwrap()
                    .iter()
                    .map(|t| {
                        let table = t.vefc.as_ref().unwrap().codon_table;
                        (t.chr.clone(), t.stable_id.to_string(), table)
                    })
                    .collect()
            };
            let pair = format!("tables {first} and {second}");
            assert_eq!(tables("1"), loaded_as("1", first), "{pair}");
            assert_eq!(tables("2"), loaded_as("2", second), "{pair}");
            assert_eq!(tables("1"), loaded_as("1", first), "{pair}");
        }
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
        let facts = *convert_transcript(jt, "1", &mut AttributeCodes::default())
            .unwrap()
            .facts();
        assert!(facts.has_frameshift_intron);
        assert_eq!(
            facts.has_frameshift_intron,
            facts.vefc_has_frameshift_intron
        );

        let without_vefc = r#"{"stable_id":"ENST2","gene_stable_id":"ENSG2","start":100,"end":300,
            "strand":-1,"biotype":"protein_coding","source":"ensembl",
            "exons":[{"start":215,"end":300},{"start":100,"end":200}]}"#;
        let jt: JsonTranscript = serde_json::from_str(without_vefc).unwrap();
        let facts = *convert_transcript(jt, "1", &mut AttributeCodes::default())
            .unwrap()
            .facts();
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
        let vefc = convert_transcript(jt, "1", &mut AttributeCodes::default())
            .unwrap()
            .vefc
            .unwrap();
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
        assert!(convert_transcript(jt, "1", &mut AttributeCodes::default())
            .unwrap()
            .vefc
            .unwrap()
            .seq_edits
            .is_empty());
    }

    /// A one-exon coding transcript on chromosome 21 at 1000-1029, ten codons
    /// `ATG GCT x8 TGA`, in the JSON the cache converter writes, with
    /// `translation_extra` spliced into its translation object.
    fn coding_transcript_json(translation_extra: &str) -> String {
        format!(
            r#"{{"stable_id":"ENST1","gene_stable_id":"ENSG1","chr":"21","start":1000,"end":1029,
            "strand":1,"biotype":"protein_coding","source":"ensembl","is_canonical":1,
            "exons":[{{"start":1000,"end":1029,"rank":1,"phase":0,"end_phase":0}}],
            "translation":{{"stable_id":"ENSP1","start":1,"end":30{translation_extra}}},
            "variation_effect_feature_cache":{{"codon_table":1,
              "translateable_seq":"ATGGCTGCTGCTGCTGCTGCTGCTGCTTGA","peptide":"MAAAAAAAA",
              "mapper":{{"start_phase":0,"cdna_coding_start":1,"cdna_coding_end":30,
                "pairs":[{{"from_start":1,"from_end":30,"to_start":1000,"to_end":1029,"ori":1}}]}}}}}}"#
        )
    }

    /// The default-format rows for a missense SNV at 21:1004 against one
    /// transcript, under `options`.
    fn missense_rows_under(
        transcript_json: &str,
        variations: Option<Arc<VariationSource>>,
        options: FieldOptions,
    ) -> String {
        let jt: JsonTranscript = serde_json::from_str(transcript_json).unwrap();
        let mut raw = HashMap::new();
        raw.insert(
            "21".to_string(),
            vec![convert_transcript(jt, "21", &mut AttributeCodes::default()).unwrap()],
        );
        let transcripts = Arc::new(LazyTranscriptIndexes::from_loaded(raw, 5_000, 5_000));
        crate::runner::default_rows_for_test(
            "21\t1004\tsnv\tC\tA\t.\t.\t.",
            transcripts,
            variations,
            options,
        )
    }

    /// The default-format rows for a missense SNV at 21:1004 against one transcript.
    fn missense_rows(transcript_json: &str, variations: Option<Arc<VariationSource>>) -> String {
        let mut options = FieldOptions::default();
        options.check_existing = variations.is_some();
        missense_rows_under(transcript_json, variations, options)
    }

    /// The fixed columns of a variation record on the region layout, which has
    /// no chromosome column.
    const REGION_COLUMNS: &str = "variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms";

    /// The metadata of a region-layout cache whose records carry the fixed
    /// columns then `extra_columns`.
    fn region_cache_info(extra_columns: &[&str]) -> CacheInfo {
        let mut info = CacheInfo::default();
        info.variation_cols = REGION_COLUMNS
            .split(',')
            .chain(extra_columns.iter().copied())
            .map(str::to_string)
            .collect();
        info
    }

    /// Writes `line` as the one record of the gzip region file at `path`.
    fn write_region_file(path: &Path, line: &str) {
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(path).unwrap(),
            flate2::Compression::default(),
        );
        writeln!(encoder, "{line}").unwrap();
        encoder.finish().unwrap();
    }

    /// A variation source over chromosome 21 whose one region file, in the
    /// cache's own `<start>-<end>_var.gz` layout, holds `line`, a record under
    /// the fixed columns then `extra_columns`. The directory must outlive the
    /// source, which parses the file on first use.
    fn region_source(
        extra_columns: &[&str],
        line: &str,
    ) -> (tempfile::TempDir, Arc<VariationSource>) {
        let dir = tempfile::tempdir().unwrap();
        let region_dir = dir.path().join("variations").join("21");
        std::fs::create_dir_all(&region_dir).unwrap();
        write_region_file(&region_dir.join("1-1000000_var.gz"), line);
        let source = VariationSource::open(
            dir.path(),
            &region_cache_info(extra_columns),
            false,
            Arc::new(CacheManifest::absent(dir.path())),
        )
        .unwrap()
        .expect("a region file is a variation source");
        (dir, Arc::new(source))
    }

    /// A variation source over chromosome 21 holding one record in the
    /// cache's own `<start>-<end>_var.gz` layout: `extra_columns` extends the
    /// fixed columns, `extra_fields` the record's cells in the same order. The
    /// directory must outlive the source, which parses the file on first use.
    fn variations_from(
        extra_columns: &[&str],
        extra_fields: &[&str],
    ) -> (tempfile::TempDir, Arc<VariationSource>) {
        let mut fields = vec![
            "rs1", "0", "0", "1004", "", "C/A", "1", "", "0", "", "", "", "", "",
        ];
        fields.extend_from_slice(extra_fields);
        region_source(extra_columns, &fields.join(" "))
    }

    /// The one record of such a source, as the reader parses it.
    fn only_record(source: &VariationSource) -> CachedVariation {
        let records = source.records_overlapping("21", 1004, 1004).unwrap();
        assert_eq!(records.len(), 1);
        records.into_iter().next().unwrap()
    }

    /// A variation record's `AA` and `EA` columns change nothing in the output
    /// under `--check_existing`, the flag that reads a co-located record:
    /// Ensembl VEP's `OutputFactory.pm` `%FREQUENCY_KEYS` has no AA or EA entry,
    /// and no release-115 or later cache carries the columns.
    #[test]
    fn legacy_esp_columns_change_no_output() {
        let (_without_dir, without) = variations_from(&["AFR"], &["A:0.01"]);
        let (_with_dir, with) =
            variations_from(&["AFR", "AA", "EA"], &["A:0.01", "A:0.2", "A:0.3"]);
        assert_eq!(only_record(&with).frequencies["AA"], "A:0.2");
        assert!(!only_record(&without).frequencies.contains_key("AA"));
        let transcript = coding_transcript_json("");
        let rows = missense_rows(&transcript, Some(with));
        assert!(
            rows.contains("missense_variant") && rows.contains("\trs1\t"),
            "{rows}"
        );
        assert_eq!(rows, missense_rows(&transcript, Some(without)));
    }

    /// A region line whose `clin_sig_allele` and `clinical_impact` cells are empty
    /// inside the line prints as Ensembl VEP 116.3 prints a region-layout cache
    /// (`21 10602054 ATA A` against rs149363218):
    /// `CLIN_SIG` is left out under the default `--clin_sig_allele 1`
    /// (OutputFactory.pm:1097, 1121, 1125 on the defined `''`), `CLIN_SIG=benign`
    /// prints under `--clin_sig_allele 0`, and `--clinvar_somatic_classification`
    /// writes `CLINVAR_SOMATIC_CLASSIFICATION=` (:1093). A line ending after
    /// `phenotype_or_disease` leaves both cells undefined, so `CLIN_SIG=benign`
    /// prints under the default and the classification key is not written.
    #[test]
    fn region_line_empty_clinvar_cells_print_as_ensembl_vep_does() {
        let rows_for = |line: &str, options: FieldOptions| {
            let (_dir, source) = region_source(&["AFR"], line);
            missense_rows_under(&coding_transcript_json(""), Some(source), options)
        };
        let defined_empty = "rs1 0 0 1004  C/A 1 benign 1      A:0.01";
        let trailing = "rs1 0 0 1004  C/A 1 benign 1";
        let default = {
            let mut options = FieldOptions::default();
            options.check_existing = true;
            options
        };
        let mut whole = default.clone();
        whole.clin_sig_allele = ClinSigAllele::WholeVariant;
        let mut somatic = default.clone();
        somatic.clinvar_somatic_classification = true;
        let rows = rows_for(defined_empty, default.clone());
        assert!(
            rows.contains("\trs1\t") && !rows.contains("CLIN_SIG"),
            "{rows}"
        );
        assert!(
            rows_for(defined_empty, whole).contains("CLIN_SIG=benign"),
            "the whole-variant term prints"
        );
        let rows = rows_for(defined_empty, somatic.clone());
        assert!(
            rows.contains(";CLINVAR_SOMATIC_CLASSIFICATION=;PHENO=1"),
            "the key is written with an empty value: {rows}"
        );
        assert!(rows_for(trailing, default).contains("CLIN_SIG=benign"));
        assert!(!rows_for(trailing, somatic).contains("CLINVAR_SOMATIC_CLASSIFICATION"));
    }

    /// A variation record's `minor_allele` and `minor_allele_freq` columns, which
    /// a pre-108 cache declares and no release-115 or later cache carries, print
    /// in no text column under `--check_existing` (`get_frequency_data` reads
    /// them for `--freq_pop 1KG_ALL` alone, BaseCacheVariation.pm:216-228) and
    /// are echoed in every JSON `colocated_variants` object under their own
    /// names, as `add_colocated_variant_info_JSON` copies the record's hash
    /// (OutputFactory/JSON.pm:415, :473): Ensembl VEP 116.3 writes both keys in
    /// every object on such a cache. The reader keeps the columns as it keeps
    /// every column outside the fixed set.
    #[test]
    fn minor_allele_columns_print_in_json_alone() {
        let (_without_dir, without) = variations_from(&[], &[]);
        let (_with_dir, with) =
            variations_from(&["minor_allele", "minor_allele_freq"], &["A", "0.12"]);
        assert_eq!(only_record(&with).frequencies["minor_allele"], "A");
        assert!(!only_record(&without)
            .frequencies
            .contains_key("minor_allele"));
        let transcript = coding_transcript_json("");
        let rows = missense_rows(&transcript, Some(Arc::clone(&with)));
        assert!(rows.contains("\trs1\t"), "{rows}");
        assert_eq!(rows, missense_rows(&transcript, Some(Arc::clone(&without))));

        let colocated = |source: &VariationSource| -> serde_json::Value {
            let mut variant =
                crate::vcf_parser::parse_vcf_line("21\t1004\tsnv\tC\tA\t.\t.\t.", true)
                    .expect("parse test VCF line")
                    .remove(0);
            crate::variation_matcher::annotate_colocated(
                &mut variant,
                source,
                &crate::variation_matcher::ColocationRules::default(),
                &crate::variation_source::ReadFailure::default(),
            );
            let mut options = FieldOptions::default();
            options.check_existing = true;
            vep_io::output::json::JsonOutputFormatter::new(options, Vec::new(), None)
                .format_variant_value(&variant)["colocated_variants"][0]
                .clone()
        };
        let echoed = colocated(&with);
        assert_eq!(echoed["id"], "rs1");
        assert_eq!(echoed["minor_allele"], "A");
        assert_eq!(echoed["minor_allele_freq"], 0.12);
        let plain = colocated(&without);
        assert_eq!(plain["id"], "rs1");
        assert!(plain.get("minor_allele").is_none());
        assert!(plain.get("minor_allele_freq").is_none());
    }

    /// A transcript's `source` names the set it came from; Ensembl VEP reads it
    /// only to resolve two transcripts sharing a stable ID under `--refseq` or
    /// `--merged` to the `ensembl` one (`AnnotationType/Transcript.pm:313-339` at
    /// release/116.3), and the SOURCE column prints the separate `_source_cache`
    /// key (`OutputFactory.pm:1592` at release/116.3), so at the default flags two
    /// caches differing only there annotate alike.
    #[test]
    fn transcript_source_changes_no_output() {
        let ensembl = coding_transcript_json("");
        let refseq = ensembl.replace(r#""source":"ensembl""#, r#""source":"RefSeq""#);
        assert_ne!(ensembl, refseq);
        let rows = missense_rows(&ensembl, None);
        assert!(rows.contains("missense_variant"), "{rows}");
        assert_eq!(rows, missense_rows(&refseq, None));
    }

    /// The translation object's `start_exon` and `end_exon` objects and its
    /// exon-relative offsets (the `Bio::EnsEMBL::Translation` fields the Storable
    /// cache keeps, `Pipeline/DumpVEP/Dumper/Core.pm:226-235` at release/116.3
    /// stripping only the translation's attributes, protein features, dates,
    /// cross-references and adaptor) are not read: the coding region comes from
    /// the exon coordinate mapper, so a cache with or without them annotates alike.
    #[test]
    fn translation_exon_fields_change_no_output() {
        let plain = coding_transcript_json("");
        let with_exons = coding_transcript_json(
            r#","start_exon":{"start":1000,"end":1029,"stable_id":"ENSE1"},"end_exon":{"start":1000,"end":1029,"stable_id":"ENSE1"}"#,
        );
        let shifted = plain.replace(r#""start":1,"end":30"#, r#""start":7,"end":12"#);
        assert_ne!(plain, shifted);
        let rows = missense_rows(&plain, None);
        assert!(rows.contains("missense_variant"), "{rows}");
        assert_eq!(rows, missense_rows(&with_exons, None));
        assert_eq!(rows, missense_rows(&shifted, None));
    }

    /// `_gene_symbol_id`, the primary identifier of the EntrezGene entry a
    /// RefSeq core's gene symbol is taken from, is set by
    /// `Database/Transcript.pm:549-567` at release/116.3 and read by no output
    /// column: SYMBOL_SOURCE and HGNC_ID print `_gene_symbol_source` and
    /// `_gene_hgnc_id`. So the converter's `gene_symbol_id` key loads as the
    /// transcript's identifier, read by `--transcript_filter` alone, and two
    /// shards differing only by it annotate alike under `--symbol`.
    #[test]
    fn gene_symbol_id_key_is_accepted_and_changes_nothing() {
        let symbol = r#","gene_symbol":"TP53","gene_symbol_source":"EntrezGene""#;
        let without = coding_transcript_json("").replace(
            r#""source":"ensembl""#,
            &format!(r#""source":"ensembl"{symbol}"#),
        );
        let with = without.replace(symbol, &format!(r#"{symbol},"gene_symbol_id":"7157""#));
        assert_ne!(without, with);
        let loaded = |body: &str| {
            let dir = cache_with_shards("vep_jc_symbol_id_", &[("21", &format!("[{body}]"))]);
            let mut transcripts = load_all_transcripts(dir.path().to_str().unwrap()).unwrap();
            let mut transcripts = transcripts.remove("21").unwrap();
            assert_eq!(transcripts.len(), 1);
            transcripts.pop().unwrap()
        };
        let (plain, identified) = (loaded(&without), loaded(&with));
        assert_eq!(plain.gene_symbol_id, None);
        assert_eq!(identified.gene_symbol_id.as_deref(), Some("7157"));
        let mut stripped = identified.clone();
        stripped.gene_symbol_id = None;
        assert_eq!(format!("{plain:?}"), format!("{stripped:?}"));
        let mut options = FieldOptions::default();
        options.symbol = true;
        let rows = missense_rows_under(&without, None, options.clone());
        assert!(
            rows.contains("SYMBOL=TP53") && rows.contains("SYMBOL_SOURCE=EntrezGene"),
            "{rows}"
        );
        assert_eq!(rows, missense_rows_under(&with, None, options));
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
        convert_transcript(jt, "21", &mut AttributeCodes::default())
            .unwrap_or_else(|e| panic!("{e}"))
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
    /// (`OutputFactory.pm:1619-1620` at release/116.3, the same rule the pick
    /// order reads at `:752-753`): `tsl1 (assigned to previous version 8)` is
    /// `1`, `tslNA` with or without a suffix is none, and a zero is none. The
    /// level is the digit run as text, since `$1` is a string: `tsl007` is
    /// `007` (116.3 prints `TSL=007`), `tsl00` is `00` (the string is true, so
    /// 116.3 prints `TSL=00`), and `tsl256` and `tsl999` keep their digits.
    #[test]
    fn tsl_leading_digits_with_assigned_suffix() {
        let tsl = |value: &str| transcript_from(&minimal(&attributes(&[("TSL", value)]))).tsl;
        assert_eq!(
            tsl("tsl1 (assigned to previous version 8)").as_deref(),
            Some("1")
        );
        assert_eq!(tsl("tsl5").as_deref(), Some("5"));
        assert_eq!(
            tsl("tsl3 (assigned to previous version 2)").as_deref(),
            Some("3")
        );
        assert_eq!(tsl("tslNA"), None);
        assert_eq!(tsl("tslNA (assigned to previous version 4)"), None);
        assert_eq!(tsl("tsl0"), None);
        // The first digit-bearing match is the only one read, as the regex's is.
        assert_eq!(tsl("tsl0 (assigned to previous version tsl3)"), None);
        assert_eq!(tsl("tslNA then tsl4").as_deref(), Some("4"));
        for level in ["007", "00", "255", "256", "999"] {
            let suffixed = format!("tsl{level} (assigned to previous version 7)");
            assert_eq!(tsl(&suffixed).as_deref(), Some(level), "{suffixed}");
        }
        assert_eq!(tsl_level("tsl2"), Some("2"));
        assert_eq!(tsl_level(""), None);

        // A cache writing the level as a top-level key, a number or a string of
        // digits, is read as that text; a zero there is none as well.
        assert_eq!(
            transcript_from(&minimal(r#","tsl":2"#)).tsl.as_deref(),
            Some("2")
        );
        assert_eq!(
            transcript_from(&minimal(r#","tsl":"007""#)).tsl.as_deref(),
            Some("007")
        );
        assert_eq!(transcript_from(&minimal(r#","tsl":0"#)).tsl, None);
        // The key stands in the attribute's place: a zero there is none even
        // beside an attribute with a level.
        let zero_beside_three = minimal(&format!(r#","tsl":0{}"#, attributes(&[("TSL", "tsl3")])));
        assert_eq!(transcript_from(&zero_beside_three).tsl, None);
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
            let txs = load_transcripts_for_chr(chr, &shards[chr], &unchecked()).unwrap();
            assert_eq!(txs.len(), 1);
            assert_eq!(&*txs[0].stable_id, id);
            assert_eq!(txs[0].chr, chr);
        }
    }

    /// The key set is Ensembl VEP's `valid_chromosomes`, every subdirectory
    /// not starting with a dot (CacheDir.pm:431): a chromosome directory with
    /// no shard is a chromosome of the cache holding no transcript, and a
    /// dot-prefixed directory is none, shard or not.
    #[test]
    fn empty_chromosome_directory_is_a_key_and_a_dot_directory_none() {
        let dir = cache_with_shards("vep_jc_valid_", &[("21", "[]")]);
        let transcripts = dir.path().join("transcripts");
        std::fs::create_dir(transcripts.join("GL000008.2")).unwrap();
        std::fs::create_dir(transcripts.join(".hidden")).unwrap();
        std::fs::write(transcripts.join(".hidden").join("1-1000000.json"), "[]").unwrap();
        let root = dir.path().to_str().unwrap();
        let shards = enumerate_transcript_shards(root).unwrap();
        let mut chromosomes: Vec<&String> = shards.keys().collect();
        chromosomes.sort();
        assert_eq!(chromosomes, ["21", "GL000008.2"]);
        assert!(shards["GL000008.2"].is_empty());
        let indexes = LazyTranscriptIndexes::new(root, 5_000, 5_000).unwrap();
        indexes.prewarm("GL000008.2").unwrap();
        assert!(indexes.contains_key("GL000008.2"));
        assert!(indexes.get("GL000008.2").unwrap().is_empty());
        assert!(!indexes.contains_key(".hidden"));
    }

    /// The whole-cache refusal counts shards, not chromosome directories: a
    /// cache whose only shard sits in a dot-prefixed directory holds none. The
    /// refusal is this engine's own: Ensembl VEP lists every non-dot directory as
    /// a valid chromosome whatever it holds (`CacheDir.pm:431`), reads a region's
    /// file only `if -e $file` (`BaseSerialized.pm:263`) and runs such a cache
    /// to the end, every record `intergenic_variant`, exit 0, as
    /// `docs/divergences/damaged-cache-shard.md` states.
    #[test]
    fn cache_with_shards_only_in_a_dot_directory_is_refused() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_dotonly_")
            .tempdir()
            .unwrap();
        let transcripts = dir.path().join("transcripts");
        std::fs::create_dir_all(transcripts.join("21")).unwrap();
        std::fs::create_dir_all(transcripts.join(".hidden")).unwrap();
        std::fs::write(transcripts.join(".hidden").join("1-1000000.json"), "[]").unwrap();
        let err = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap_err();
        assert!(
            CliError::from_run(&err).starts_with("ERROR: [cache] no transcript shards under"),
            "{err}"
        );
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

    /// `gene_stable_id` (VEP's `_gene_stable_id`) is carried onto the loaded
    /// transcript. A transcript whose gene has no stable id, written without
    /// the key or with a null, loads with an empty one, the value the writers
    /// print as no gene; the gene's `_gene` object is not in the cache and
    /// nothing reads it.
    #[test]
    fn gene_stable_id_loaded_and_empty_when_absent() {
        assert_eq!(
            &*transcript_from(&minimal("")).gene_stable_id,
            "ENSG00000215455"
        );
        let without = r#"{"stable_id":"NM_1","start":1,"end":9,"strand":1,"biotype":"cDNA","source":"ensembl"}"#;
        assert_eq!(&*transcript_from(without).gene_stable_id, "");
        let null = r#"{"stable_id":"NM_1","gene_stable_id":null,"start":1,"end":9,"strand":1,"biotype":"cDNA","source":"ensembl"}"#;
        assert_eq!(&*transcript_from(null).gene_stable_id, "");
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

    /// A translation the cache holds without a stable id reaches the shard
    /// without the `stable_id` key, the converter writing defined values alone.
    /// It loads with an empty id and its database id, and the transcript is
    /// otherwise the one a translation with an id gives.
    #[test]
    fn a_translation_without_a_stable_id_loads_with_an_empty_id() {
        let with_id = coding_transcript_json(r#","dbID":7001"#);
        let without_id = with_id.replacen(r#""stable_id":"ENSP1","#, "", 1);
        assert_ne!(without_id, with_id);

        let loaded = transcript_from(&without_id);
        let translation = loaded.translation.as_ref().unwrap();
        assert_eq!(translation.stable_id, "");
        assert_eq!(translation.db_id, Some(7001));

        let mut expected = transcript_from(&with_id);
        expected.translation.as_mut().unwrap().stable_id.clear();
        assert_eq!(
            serde_json::to_value(&loaded).unwrap(),
            serde_json::to_value(&expected).unwrap()
        );
    }

    /// Ensembl VEP's cache stores `coding_region_start` and `coding_region_end`
    /// undefined and computes both at run time from the translation's exons, so
    /// the loader derives them from the mapper's coding cDNA span, in genomic
    /// order on either strand, and a transcript with the two keys loads to the
    /// same `Transcript` as one without them.
    #[test]
    fn coding_region_bounds_recomputed_from_the_mapper() {
        let forward = coding_transcript_json("");
        let reverse = minimal(
            r#","translation":{"stable_id":"ENSP1","start":4,"end":12},"exons":[{"start":44538981,"end":44538995}],
            "variation_effect_feature_cache":{"translateable_seq":"ATGCCCTAA",
            "mapper":{"start_phase":-1,"cdna_coding_start":4,"cdna_coding_end":12,"pairs":[{"from_start":1,"from_end":15,"to_start":44538981,"to_end":44538995,"ori":-1}]}}"#,
        );
        // Forward: cDNA 1..30 is genomic 1000..1029. Reverse: cDNA 4 and 12 are
        // genomic 44538995 - 3 and 44538995 - 11, swapped into start-end order.
        for (json, (start, end)) in [(forward, (1000, 1029)), (reverse, (44538984, 44538992))] {
            let derived = transcript_from(&json);
            assert_eq!(
                (derived.coding_region_start, derived.coding_region_end),
                (Some(start), Some(end))
            );
            assert_eq!(
                (derived.translation_start, derived.translation_end),
                (Some(start), Some(end))
            );
            let stated = json.replacen(
                r#"{"stable_id":"#,
                r#"{"coding_region_start":1,"coding_region_end":2,"stable_id":"#,
                1,
            );
            assert_ne!(stated, json);
            assert_eq!(
                serde_json::to_value(transcript_from(&stated)).unwrap(),
                serde_json::to_value(&derived).unwrap()
            );
        }
    }

    /// A coding bound in bases the model carries and the genome lacks (a block
    /// an `_rna_edit` inserts, which no pair holds) is read at the nearest
    /// coding base a pair does hold, on either strand: the first one at or
    /// after the coding start, the last one at or before the coding end.
    /// Ensembl VEP takes both bounds from the translation's exons (ensembl
    /// core `Transcript.pm:1087-1154`), so they exist whatever the pairs
    /// hold. A coding sequence with no base on the genome has no bounds.
    #[test]
    fn a_coding_bound_no_pair_holds_is_read_at_the_nearest_paired_coding_base() {
        let bounds = |strand: i8, coding: (u64, u64), pairs: &[(u64, u64, u64, u64)]| {
            let pairs: Vec<String> = pairs
                .iter()
                .map(|(from_start, from_end, to_start, to_end)| {
                    format!(
                        r#"{{"from_start":{from_start},"from_end":{from_end},
                        "to_start":{to_start},"to_end":{to_end},"ori":{strand}}}"#
                    )
                })
                .collect();
            let transcript = transcript_from(&format!(
                r#"{{"stable_id":"NM_1","gene_stable_id":"G1","start":1000,"end":1029,
                "strand":{strand},"biotype":"protein_coding","source":"refseq",
                "translation":{{"stable_id":"NP_1","start":1,"end":30}},
                "exons":[{{"start":1000,"end":1029}}],
                "variation_effect_feature_cache":{{"mapper":{{"start_phase":0,
                "cdna_coding_start":{},"cdna_coding_end":{},"pairs":[{}]}}}}}}"#,
                coding.0,
                coding.1,
                pairs.join(",")
            ));
            assert_eq!(
                (transcript.translation_start, transcript.translation_end),
                (transcript.coding_region_start, transcript.coding_region_end)
            );
            (transcript.coding_region_start, transcript.coding_region_end)
        };
        let whole = Some(1000).zip(Some(1029));
        let both = |bounds: (Option<u64>, Option<u64>)| bounds.0.zip(bounds.1);

        // Six bases before the first genomic base, the coding start among them.
        assert_eq!(both(bounds(1, (1, 36), &[(7, 36, 1000, 1029)])), whole);
        assert_eq!(both(bounds(-1, (1, 36), &[(7, 36, 1000, 1029)])), whole);
        // Six bases after the last genomic base, the coding end among them.
        assert_eq!(both(bounds(1, (1, 36), &[(1, 30, 1000, 1029)])), whole);
        assert_eq!(
            both(bounds(-1, (4, 36), &[(1, 30, 1000, 1029)])),
            Some((1000, 1026))
        );
        // Six bases between two pairs, the coding start or the coding end
        // among them.
        let split = [(1, 10, 1000, 1009), (17, 36, 1010, 1029)];
        assert_eq!(both(bounds(1, (12, 36), &split)), Some((1010, 1029)));
        assert_eq!(both(bounds(1, (1, 13), &split)), Some((1000, 1009)));
        // Both bounds in blocks, and both on the genome.
        assert_eq!(both(bounds(1, (1, 42), &[(7, 36, 1000, 1029)])), whole);
        assert_eq!(
            both(bounds(1, (9, 34), &[(7, 36, 1000, 1029)])),
            Some((1002, 1027))
        );
        // No coding base on the genome.
        assert_eq!(bounds(1, (1, 6), &[(7, 36, 1000, 1029)]), (None, None));
        assert_eq!(bounds(1, (37, 42), &[(7, 36, 1000, 1029)]), (None, None));
        // A coding end before the coding start, each in its own pair: the
        // span between the two.
        let thirds = [
            (1, 10, 1000, 1009),
            (11, 20, 1010, 1019),
            (21, 30, 1020, 1029),
        ];
        assert_eq!(both(bounds(1, (25, 5), &thirds)), Some((1004, 1024)));
        // A pair whose end precedes its start holds no base.
        let empty = [(1, 9, 1000, 1008), (20, 10, 1009, 1019)];
        assert_eq!(both(bounds(1, (5, 30), &empty)), Some((1004, 1008)));
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

    /// Ensembl VEP's cache dumper deletes a transcript's `description` when it
    /// is defined, one of the keys it drops with `display_xref`, `external_name`
    /// and `status` (`Pipeline/DumpVEP/Dumper/Core.pm:155-173`), so no cache
    /// carries a value under the key: a record that does loads with the key
    /// ignored and the transcript carries no such field.
    #[test]
    fn description_key_is_ignored() {
        let tx = transcript_from(&minimal(
            r#","description":"keratin associated protein 10-1 [Source:HGNC Symbol;Acc:HGNC:22966]""#,
        ));
        assert_eq!(&*tx.stable_id, "ENST00000400375");
        let fields = serde_json::to_value(&tx).unwrap();
        assert!(fields.get("description").is_none(), "{fields}");
    }

    /// On a parse keeping exon ids they travel on the transcript in rank
    /// order, whatever order the shard lists the exons in; an exon lacking an
    /// id holds the empty id, and a transcript whose exons all lack one, or
    /// with no exons, carries none; the ids are the loaded exons' companions,
    /// so each shares its exon's index.
    #[test]
    fn exon_ids_are_kept_in_rank_order_with_the_empty_id_for_an_exon_without_one() {
        let exon = |rank: u32, id: Option<&str>| JsonExon {
            stable_id: id.map(Arc::from),
            start: 100 * i64::from(rank),
            end: 100 * i64::from(rank) + 50,
            rank: Some(rank),
            phase: -1,
            end_phase: -1,
        };
        let ids = exon_ids_in_rank_order(&[
            exon(2, Some("ENSE2")),
            exon(1, Some("ENSE1")),
            exon(3, Some("ENSE3")),
        ])
        .unwrap();
        assert_eq!(
            ids.iter().map(|id| &**id).collect::<Vec<_>>(),
            ["ENSE1", "ENSE2", "ENSE3"]
        );
        let some = exon_ids_in_rank_order(&[exon(2, None), exon(1, Some("ENSE1"))]).unwrap();
        assert_eq!(
            some.iter().map(|id| &**id).collect::<Vec<_>>(),
            ["ENSE1", ""]
        );
        assert!(exon_ids_in_rank_order(&[exon(1, None), exon(2, None)]).is_none());
        assert!(exon_ids_in_rank_order(&[]).is_none());
    }

    /// A shard parse keeping exon ids puts them on the transcript, and the
    /// option ends with the parse: the next parse on the thread, asked for
    /// none, loads the same shard without them.
    #[test]
    fn a_parse_keeps_exon_ids_for_its_own_duration_alone() {
        let shard = format!(
            "[{}]",
            minimal(
                r#","exons":[{"stable_id":"ENSE00001472863","start":44540001,"end":44540195,"rank":2,"phase":-1,"end_phase":-1},{"stable_id":"ENSE00001472862","start":44538981,"end":44539100,"rank":1,"phase":-1,"end_phase":-1}]"#
            )
        );
        let exon_ids = |keep: bool| {
            let mut parsed = parse_shard_transcripts(shard.as_bytes(), keep, true).unwrap();
            let jt = parsed.pop().unwrap();
            convert_transcript(jt, "21", &mut AttributeCodes::default())
                .unwrap()
                .exon_ids
        };
        let kept = exon_ids(true).expect("the ids are kept");
        assert_eq!(
            kept.iter().map(|id| &**id).collect::<Vec<_>>(),
            ["ENSE00001472862", "ENSE00001472863"]
        );
        assert!(!KEEP_EXON_IDS.get(), "the option ends with the parse");
        assert!(exon_ids(false).is_none());
    }

    /// A shard parse without predictions drops the SIFT and PolyPhen matrices,
    /// and the option ends with the parse: the next parse on the thread keeps
    /// them, as every load does unless asked otherwise.
    #[test]
    fn a_parse_drops_predictions_for_its_own_duration_alone() {
        let shard = format!(
            "[{}]",
            minimal(
                r#","variation_effect_feature_cache":{"protein_function_predictions":{"sift":{"analysis":"sift","predictions_data":"AAAA"}}}"#
            )
        );
        let predictions = |keep: bool| {
            let mut parsed = parse_shard_transcripts(shard.as_bytes(), false, keep).unwrap();
            let jt = parsed.pop().unwrap();
            convert_transcript(jt, "21", &mut AttributeCodes::default())
                .unwrap()
                .vefc
                .unwrap()
                .protein_function_predictions
        };
        assert!(predictions(false).is_none());
        assert!(KEEP_PREDICTIONS.get(), "the option ends with the parse");
        assert!(predictions(true).is_some_and(|p| p.sift.is_some()));
    }

    /// The keys no output or consequence reads (the transcript's `description`,
    /// each exon's `stable_id` on a parse not keeping exon ids) are neither
    /// loaded nor required: a record carrying them and the same record without
    /// them load to equal transcripts, and the loaded exon carries no id.
    #[test]
    fn dropped_fields_are_ignored_not_required() {
        let exons_with_ids = r#","exons":[{"stable_id":"ENSE00001472862","start":44538981,"end":44539100,"rank":1,"phase":-1,"end_phase":-1},{"stable_id":"ENSE00001472863","start":44540001,"end":44540195,"rank":2,"phase":-1,"end_phase":-1}]"#;
        let exons_without_ids = r#","exons":[{"start":44538981,"end":44539100,"rank":1,"phase":-1,"end_phase":-1},{"start":44540001,"end":44540195,"rank":2,"phase":-1,"end_phase":-1}]"#;
        let with = transcript_from(&minimal(&format!(
            r#","description":"keratin associated protein 10-1"{exons_with_ids}"#
        )));
        let without = transcript_from(&minimal(exons_without_ids));
        assert_eq!(with.exons.len(), 2);
        assert_eq!(
            serde_json::to_value(&with).unwrap(),
            serde_json::to_value(&without).unwrap()
        );
        let exon = serde_json::to_value(&with.exons[0]).unwrap();
        assert!(exon.get("stable_id").is_none(), "{exon}");
        assert_eq!(exon["start"], 44538981);
    }

    /// Two attributes of one shard carrying the same code share one
    /// allocation, across transcripts as within one, while different codes
    /// and the attribute values stay their own.
    #[test]
    fn attribute_codes_are_interned_once_per_load() {
        let first = minimal(&attributes(&[
            ("TSL", "tsl1"),
            ("gencode_basic", "GENCODE basic"),
        ]));
        let second = minimal(&attributes(&[("TSL", "tsl5"), ("appris", "principal1")]));
        let shard = format!(r#"{{"21":[{first},null,{second}]}}"#);
        let mut codes = AttributeCodes::default();
        let txs: Vec<Transcript> = parse_shard_transcripts(shard.as_bytes(), false, true)
            .unwrap()
            .into_iter()
            .map(|jt| convert_transcript(jt, "21", &mut codes).unwrap())
            .collect();
        assert_eq!(txs.len(), 2);
        let (a, b) = (&txs[0].attributes, &txs[1].attributes);
        assert_eq!((&*a[0].code, &*b[0].code), ("TSL", "TSL"));
        assert!(Arc::ptr_eq(&a[0].code, &b[0].code));
        assert!(!Arc::ptr_eq(&a[0].code, &a[1].code));
        assert_eq!((a[0].value.as_str(), b[0].value.as_str()), ("tsl1", "tsl5"));
        assert_eq!(codes.0.len(), 3);
        assert_eq!(
            (txs[0].tsl.as_deref(), txs[1].tsl.as_deref()),
            (Some("1"), Some("5"))
        );
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
        assert_eq!((hit.label(), hit.score), ("deleterious", 0.01));
        assert!(sift.lookup(1, b'A', AnalysisType::Sift).is_none());
        assert!(predictions.polyphen_humvar.is_none());
        assert!(predictions.polyphen_humdiv.is_none());
    }

    /// A `predictions_data` that is empty or not base64 loads as no matrix,
    /// as the key's absence does; base64 of bytes that are no gzip stream
    /// loads as a matrix that stays as stored until its first lookup and
    /// then predicts nothing at any residue.
    #[test]
    fn malformed_predictions_data_predicts_nothing() {
        use vep_core::prediction::AnalysisType;

        let sift_of = |field: &str| {
            transcript_from(&minimal(&format!(
                r#","variation_effect_feature_cache":{{"protein_function_predictions":{{"sift":{{"analysis":"sift"{field}}}}}}}"#
            )))
            .vefc
            .unwrap()
            .protein_function_predictions
            .unwrap()
            .sift
        };
        assert!(sift_of("").is_none());
        assert!(sift_of(r#","predictions_data":"""#).is_none());
        assert!(sift_of(r#","predictions_data":"not base64!""#).is_none());
        // `VkVQ` is the base64 of the bare `VEP` header, no gzip stream around it.
        let stored = sift_of(r#","predictions_data":"VkVQ""#).unwrap();
        assert!(!stored.is_inflated());
        assert!(stored.lookup(1, b'A', AnalysisType::Sift).is_none());
        assert!(stored.is_inflated() && stored.predictions_data().is_none());
    }

    /// A chromosome directory reached through a symlink loads like a real one
    /// in both walks: a skipped link would leave the transcript cache refused
    /// as shard-less and the variations silently absent.
    #[cfg(unix)]
    #[test]
    fn symlinked_chromosome_directory_is_walked() {
        let dir = tempfile::Builder::new()
            .prefix("vep_jc_symlink_")
            .tempdir()
            .unwrap();
        let real = dir.path().join("real_21");
        std::fs::create_dir_all(&real).unwrap();
        std::fs::write(real.join("1-1000000.json"), b"[]").unwrap();
        let transcripts = dir.path().join("transcripts");
        std::fs::create_dir_all(&transcripts).unwrap();
        std::os::unix::fs::symlink(&real, transcripts.join("21")).unwrap();
        let shards = enumerate_transcript_shards(dir.path().to_str().unwrap()).unwrap();
        assert_eq!(
            shards["21"],
            vec![transcripts.join("21").join("1-1000000.json")]
        );
        let real_variations = dir.path().join("real_v21");
        std::fs::create_dir_all(&real_variations).unwrap();
        write_region_file(
            &real_variations.join("1-1000000_var.gz"),
            "rs123 0 0 100  A/G 1  0",
        );
        let variations = dir.path().join("variations");
        std::fs::create_dir_all(&variations).unwrap();
        std::os::unix::fs::symlink(&real_variations, variations.join("21")).unwrap();
        let source = VariationSource::open(
            dir.path(),
            &region_cache_info(&[]),
            false,
            Arc::new(CacheManifest::absent(dir.path())),
        )
        .unwrap()
        .expect("the linked chromosome directory is a variation source");
        let records = source.records_overlapping("21", 100, 100).unwrap();
        assert_eq!(records[0].variation_name, "rs123");
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

    /// A coding transcript at 21:1000 whose one exon ends `held` bases into
    /// the stop codon `stop`, as a cache holds it: four sense codons, an
    /// `_rna_edit` replacing the exon's last `held` bases by the codon, and
    /// the two pairs Ensembl's mapper builds for the edit, the second placing
    /// the replaced genomic bases on the codon's last bases. `strand` is 1 or
    /// -1; the codon table is the vertebrate mitochondrial one.
    fn completed_stop_transcript(held: u64, stop: &str, strand: i8) -> Transcript {
        transcript_from(&completed_stop_record(held, stop, strand))
    }

    /// The cache record [`completed_stop_transcript`] reads.
    fn completed_stop_record(held: u64, stop: &str, strand: i8) -> String {
        let coding = format!("ATG{}{stop}", "GCT".repeat(3));
        let end = 1011 + held;
        let (before, replaced) = if strand == 1 {
            ((1000, 1011), (1012, end))
        } else {
            ((1000 + held, end), (1000, 999 + held))
        };
        format!(
            r#"{{"stable_id":"NM_1","gene_stable_id":"G1","start":1000,"end":{end},
            "strand":{strand},"biotype":"protein_coding","source":"refseq",
            "attributes":[{{"code":"_rna_edit","value":"13 {edit_end} {stop}"}}],
            "exons":[{{"start":1000,"end":{end}}}],
            "translation":{{"stable_id":"NP_1","start":1,"end":{edit_end}}},
            "variation_effect_feature_cache":{{"codon_table":2,"translateable_seq":"{coding}",
            "mapper":{{"start_phase":0,"cdna_coding_start":1,"cdna_coding_end":15,"pairs":[
            {{"from_start":1,"from_end":12,"to_start":{},"to_end":{},"ori":{strand}}},
            {{"from_start":{},"from_end":15,"to_start":{},"to_end":{},"ori":{strand}}}]}}}}}}"#,
            before.0,
            before.1,
            16 - held,
            replaced.0,
            replaced.1,
            edit_end = 12 + held,
        )
    }

    /// The mapper pairs of `transcript` as `(cDNA start, cDNA end, genomic
    /// start, genomic end, orientation)`.
    fn pair_tuples(transcript: &Transcript) -> Vec<(u64, u64, u64, u64, i8)> {
        let mapper = transcript.vefc.as_ref().unwrap().mapper.as_ref().unwrap();
        mapper
            .exon_coord_mapper
            .pairs
            .iter()
            .map(|p| (p.from_start, p.from_end, p.to_start, p.to_end, p.ori))
            .collect()
    }

    /// Where the reference sequence shows the replaced genomic bases to be the
    /// stop codon's first and no other placement of them to be an identity, the
    /// last pair moves onto the codon's first bases and joins the pair before
    /// it, on either strand; the coding end and the genomic coding bounds stay
    /// as loaded, and the codon's remaining bases have no genomic partner.
    #[test]
    fn a_completed_stop_codon_pairs_with_the_bases_the_genome_holds() {
        for (held, strand, genome) in [(1, 1, "T"), (2, 1, "TA"), (1, -1, "A"), (2, -1, "TA")] {
            let mut transcript = completed_stop_transcript(held, "TAA", strand);
            let replaced = if strand == 1 {
                (1012, 1011 + held)
            } else {
                (1000, 999 + held)
            };
            let bounds = (transcript.coding_region_start, transcript.coding_region_end);
            assert_eq!(bounds, (Some(1000), Some(1011 + held)));
            let changed = pair_completed_stop_codon(&mut transcript, |chr, start, end| {
                assert_eq!((chr, start, end), ("21", replaced.0, replaced.1));
                Some(genome.as_bytes().to_vec())
            });
            assert!(changed, "{held} {strand}");
            assert_eq!(
                pair_tuples(&transcript),
                [(1, 12 + held, 1000, 1011 + held, strand)],
                "{held} {strand}"
            );
            let mapper = transcript.vefc.as_ref().unwrap().mapper.as_ref().unwrap();
            assert_eq!(mapper.cdna_coding_end, 15);
            assert_eq!(mapper.exon_coord_mapper.by_genomic_start(), [0]);
            assert_eq!(
                (transcript.coding_region_start, transcript.coding_region_end),
                bounds
            );
        }
    }

    /// A set that was indexed before it met the reference sequence is paired
    /// in a copy: the copy holds the moved pair, shares the transcript it did
    /// not move and reads a shard as its index reads it. The set it was made
    /// from keeps the cache's pairs, and a paired set moves nothing.
    #[test]
    fn a_set_indexed_without_the_reference_is_paired_in_a_copy() {
        let shard = format!("[{},{}]", completed_stop_record(1, "TAA", 1), minimal(""));
        let dir = cache_with_shards("vep_jc_paired_", &[("21", &shard)]);
        let mut bases = vec![b'C'; 1100];
        bases[999..1012].copy_from_slice(b"ATGGCTGCTGCTT");
        let fasta = dir.path().join("reference.fa");
        std::fs::write(&fasta, [b">21\n".as_slice(), &bases, b"\n"].concat()).unwrap();
        let index = dir.path().join("reference.fa.fai");
        std::fs::write(index, "21\t1100\t4\t1100\t1101\n").unwrap();
        let reference = Arc::new(vep_fasta::IndexedFasta::from_path(&fasta).unwrap());

        let set = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
        set.prewarm("21").unwrap();
        let paired = set.paired_with(&reference).expect("a pair moves");
        let given = set.get("21").unwrap().transcripts();
        let copy = paired.get("21").unwrap().transcripts();
        let completed = given.iter().position(|t| &*t.stable_id == "NM_1").unwrap();
        let other = 1 - completed;
        let whole = [(1, 13, 1000, 1012, 1)];
        assert_eq!(pair_tuples(&copy[completed]), whole);
        assert_eq!(pair_tuples(&given[completed]).len(), 2);
        assert!(Arc::ptr_eq(&given[other], &copy[other]));
        assert!(paired.paired_with(&reference).is_none());

        let shard = &paired.shards_of("21").unwrap()[0];
        let read = |set: &LazyTranscriptIndexes| {
            let transcripts = set.load_shard("21", shard).unwrap();
            let completed = transcripts.iter().find(|t| &*t.stable_id == "NM_1");
            pair_tuples(completed.unwrap())
        };
        assert_eq!(read(&paired), whole);
        assert_eq!(read(&set).len(), 2);
    }

    /// An edit's span counts the unedited cDNA, so an edit before the
    /// completing one moves it by the bases it adds: with three bases inserted
    /// after cDNA 3 the edit at unedited cDNA 13 completes the codon at cDNA
    /// 16, and its pair joins the pair before it.
    #[test]
    fn a_completed_stop_codon_after_an_earlier_edit_is_paired_alike() {
        let mut transcript = completed_stop_transcript(1, "TAA", 1);
        transcript.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "4 3 GCT".into(),
        });
        let vefc = transcript.vefc.as_mut().unwrap();
        vefc.translateable_seq = Some(format!("ATG{}TAA", "GCT".repeat(4)));
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 18;
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: 1,
        };
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 3, 1000, 1002),
            pair(7, 15, 1003, 1011),
            pair(18, 18, 1012, 1012),
        ]);
        assert!(pair_completed_stop_codon(&mut transcript, |_, _, _| Some(
            b"T".to_vec()
        )));
        assert_eq!(
            pair_tuples(&transcript),
            [(1, 3, 1000, 1002, 1), (7, 16, 1003, 1012, 1)]
        );
    }

    /// A placement the reference sequence does not show keeps the cache's
    /// pairs: the genome's base is not the codon's first; it is the codon's
    /// first and its last, two identities (`AGA` is a stop codon of the
    /// vertebrate mitochondrial code); the reference has no such span, or
    /// fewer bases than the edit replaces.
    #[test]
    fn a_stop_codon_the_reference_does_not_place_keeps_the_cache_pairs() {
        let unchanged = |mut transcript: Transcript, genome: Option<&str>| {
            let before = pair_tuples(&transcript);
            let changed = pair_completed_stop_codon(&mut transcript, |_, _, _| {
                genome.map(|bases| bases.as_bytes().to_vec())
            });
            assert!(!changed, "{genome:?}");
            assert_eq!(pair_tuples(&transcript), before, "{genome:?}");
        };
        unchanged(completed_stop_transcript(1, "TAA", 1), Some("A"));
        unchanged(completed_stop_transcript(1, "TAA", 1), Some("G"));
        unchanged(completed_stop_transcript(2, "TAA", 1), Some("AA"));
        unchanged(completed_stop_transcript(1, "AGA", 1), Some("A"));
        unchanged(completed_stop_transcript(1, "TAA", 1), None);
        unchanged(completed_stop_transcript(2, "TAA", 1), Some("T"));
    }

    /// A transcript the cache does not show to end inside a stop codon an edit
    /// completes keeps its pairs, and the reference sequence is not read for
    /// it: the edit inserts before the exon's last base and replaces nothing;
    /// the codon it completes is no stop codon; an earlier edit moves its span
    /// off the final codon; a 3' UTR follows the codon; the transcript runs on
    /// past the replaced bases; the transcript has no edit at all. A coding end
    /// and a start phase at the ends of their ranges are read without a failure.
    #[test]
    fn a_transcript_of_another_shape_is_not_read_as_a_completed_stop_codon() {
        let not_read = |transcript: &mut Transcript, shape: &str| {
            let before = pair_tuples(transcript);
            let changed = pair_completed_stop_codon(transcript, |_, _, _| {
                panic!("{shape}: the reference sequence is read")
            });
            assert!(!changed, "{shape}");
            assert_eq!(pair_tuples(transcript), before, "{shape}");
        };
        let edited = || completed_stop_transcript(1, "TAA", 1);

        let mut inserting = edited();
        inserting.attributes[0].value = "13 12 TA".into();
        not_read(&mut inserting, "an inserting edit");

        let mut sense = edited();
        sense.attributes[0].value = "13 13 GCT".into();
        sense.vefc.as_mut().unwrap().translateable_seq = Some("ATG".to_string() + &"GCT".repeat(4));
        not_read(&mut sense, "a sense codon");

        let mut moved = edited();
        moved.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "4 3 GCT".into(),
        });
        not_read(&mut moved, "an earlier edit");

        let mut with_utr = edited();
        with_utr
            .vefc
            .as_mut()
            .unwrap()
            .mapper
            .as_mut()
            .unwrap()
            .cdna_coding_end = 12;
        not_read(&mut with_utr, "a 3' UTR");

        let mut longer = edited();
        longer.end += 100;
        not_read(&mut longer, "a longer transcript");

        let mut unedited = transcript_from(&coding_transcript_json(""));
        not_read(&mut unedited, "no edit");

        let mut at_range_ends = edited();
        let mapper = at_range_ends
            .vefc
            .as_mut()
            .unwrap()
            .mapper
            .as_mut()
            .unwrap();
        mapper.start_phase = i8::MAX;
        mapper.cdna_coding_end = u64::MAX;
        not_read(
            &mut at_range_ends,
            "a coding end and a phase at their range ends",
        );
    }

    /// A coding transcript as the converter writes a dumped one: the stored
    /// span, the stored exons with the introns between them, and the mapper
    /// pairs as `(cDNA start, cDNA end, genomic start, genomic end)`.
    fn stored_transcript(
        start: i64,
        end: u64,
        ori: i8,
        exons: &[(i64, i64)],
        pairs: &[(u64, u64, u64, u64)],
    ) -> Result<Transcript, String> {
        let exon_list: Vec<serde_json::Value> = exons
            .iter()
            .map(|(start, end)| serde_json::json!({"start": start, "end": end}))
            .collect();
        let introns: Vec<serde_json::Value> = exons
            .windows(2)
            .map(|pair| serde_json::json!({"start": pair[0].1 + 1, "end": pair[1].0 - 1}))
            .collect();
        let cdna_end = pairs.iter().map(|pair| pair.1).max().unwrap();
        let pair_list: Vec<serde_json::Value> = pairs
            .iter()
            .map(|(from_start, from_end, to_start, to_end)| {
                serde_json::json!({
                    "from_start": from_start,
                    "from_end": from_end,
                    "to_start": to_start,
                    "to_end": to_end,
                    "ori": ori,
                })
            })
            .collect();
        let record = serde_json::json!({
            "stable_id": "T1",
            "gene_stable_id": "G1",
            "start": start,
            "end": end,
            "strand": ori,
            "biotype": "protein_coding",
            "source": "ensembl",
            "exons": exon_list.clone(),
            "variation_effect_feature_cache": {
                "translateable_seq": "ATG".repeat((cdna_end / 3) as usize),
                "introns": introns,
                "sorted_exons": exon_list,
                "mapper": {
                    "start_phase": 0,
                    "cdna_coding_start": 1,
                    "cdna_coding_end": cdna_end,
                    "pairs": pair_list,
                },
            },
        });
        let jt: JsonTranscript = serde_json::from_value(record).unwrap();
        convert_transcript(jt, "ring", &mut AttributeCodes::default())
    }

    /// The transcript's span, exons, introns and mapper pairs' genomic spans,
    /// its coding region and the length of its circular sequence.
    #[allow(clippy::type_complexity)]
    fn placement(
        transcript: &Transcript,
    ) -> (
        (u64, u64),
        Vec<(u64, u64)>,
        Vec<(u64, u64)>,
        Vec<(u64, u64)>,
        (Option<u64>, Option<u64>),
        Option<u64>,
    ) {
        let vefc = transcript.vefc.as_ref().unwrap();
        (
            (transcript.start, transcript.end),
            transcript.exons.iter().map(|e| (e.start, e.end)).collect(),
            transcript
                .introns
                .iter()
                .map(|i| (i.start, i.end))
                .collect(),
            vefc.mapper
                .as_ref()
                .unwrap()
                .exon_coord_mapper
                .pairs
                .iter()
                .map(|p| (p.to_start, p.to_end))
                .collect(),
            (transcript.coding_region_start, transcript.coding_region_end),
            transcript.circular_length,
        )
    }

    /// A gene of 150 bases from 11911 across the origin of a sequence of
    /// 12,000 bases to base 60 is dumped with the start -89 on a pair that
    /// keeps 11911 and 60 (ensembl core `DBSQL/TranscriptAdaptor.pm:1941-1947`,
    /// `TranscriptMapper.pm:160-161`). It loads at 11911-12060 on either
    /// strand, with the sequence's length.
    #[test]
    fn a_transcript_stored_across_the_origin_runs_on_past_the_last_base() {
        for ori in [1, -1] {
            let transcript =
                stored_transcript(-89, 60, ori, &[(-89, 60)], &[(1, 150, 11_911, 60)]).unwrap();
            assert_eq!(
                placement(&transcript),
                (
                    (11_911, 12_060),
                    vec![(11_911, 12_060)],
                    vec![],
                    vec![(11_911, 12_060)],
                    (Some(11_911), Some(12_060)),
                    Some(12_000),
                ),
                "strand {ori}"
            );
        }
    }

    /// A dump can hold the exon of such a transcript past the last base
    /// already, under a transcript stored below base 1; the exon stays where
    /// it is.
    #[test]
    fn an_exon_stored_past_the_last_base_stays_where_it_is() {
        let transcript =
            stored_transcript(-89, 60, -1, &[(11_911, 12_060)], &[(1, 150, 11_911, 60)]).unwrap();
        assert_eq!(
            placement(&transcript),
            (
                (11_911, 12_060),
                vec![(11_911, 12_060)],
                vec![],
                vec![(11_911, 12_060)],
                (Some(11_911), Some(12_060)),
                Some(12_000),
            )
        );
    }

    /// Two exons either side of the origin, the intron across it: the exon
    /// among the last bases is stored below base 1 and its pair keeps its
    /// database span, the exon among the first bases is stored where it lies.
    /// Both load in one run past the last base.
    #[test]
    fn exons_on_both_sides_of_the_origin_load_in_one_run() {
        let transcript = stored_transcript(
            -89,
            60,
            1,
            &[(-89, -40), (11, 60)],
            &[(1, 50, 11_911, 11_960), (51, 100, 11, 60)],
        )
        .unwrap();
        assert_eq!(
            placement(&transcript),
            (
                (11_911, 12_060),
                vec![(11_911, 11_960), (12_011, 12_060)],
                vec![(11_961, 12_010)],
                vec![(11_911, 11_960), (12_011, 12_060)],
                (Some(11_911), Some(12_060)),
                Some(12_000),
            )
        );
    }

    /// A start below base 1 on a transcript none of whose pairs begins past
    /// its stored end names no sequence length, and the record is refused.
    #[test]
    fn a_start_below_one_without_a_pair_across_the_origin_is_refused() {
        let refused = stored_transcript(-89, 60, 1, &[(-89, 60)], &[(1, 150, 1, 60)]).unwrap_err();
        assert_eq!(
            refused,
            "transcript T1: its start -89 lies below base 1 and no mapper pair places it \
             across the origin of a circular sequence"
        );
    }

    /// A pair across the origin that starts at the largest coordinate a record
    /// can hold names a sequence of that length, and its end, moved one length
    /// on, stays at that coordinate.
    #[test]
    fn a_pair_across_the_origin_at_the_largest_coordinate_is_placed_without_overflow() {
        let record = serde_json::json!({
            "stable_id": "T1",
            "gene_stable_id": "G1",
            "start": 0,
            "end": 10,
            "strand": 1,
            "biotype": "protein_coding",
            "source": "ensembl",
            "variation_effect_feature_cache": {
                "mapper": {
                    "pairs": [{
                        "from_start": 1,
                        "from_end": 6,
                        "to_start": u64::MAX,
                        "to_end": 5,
                        "ori": 1,
                    }],
                },
            },
        });
        let mut jt: JsonTranscript = serde_json::from_value(record).unwrap();
        assert_eq!(place_across_origin(&mut jt), Ok(Some(u64::MAX)));
        let mapper = jt.vefc.unwrap().mapper.unwrap();
        let pair = &mapper.pairs[0];
        assert_eq!((pair.to_start, pair.to_end), (u64::MAX, u64::MAX));
    }

    /// A transcript that starts at base 1 is on no origin and loads as stored.
    #[test]
    fn a_transcript_that_starts_at_base_one_loads_as_stored() {
        let transcript = stored_transcript(1, 150, 1, &[(1, 150)], &[(1, 150, 1, 150)]).unwrap();
        assert_eq!(
            placement(&transcript),
            (
                (1, 150),
                vec![(1, 150)],
                vec![],
                vec![(1, 150)],
                (Some(1), Some(150)),
                None,
            )
        );
    }
}
