// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a `--custom` source leaves on a variant: the records it matched and the
//! summary statistics over them, as Ensembl VEP keeps them on the variation
//! feature (`_custom_annotations` and `_custom_annotations_stats`,
//! `AnnotationSource/File.pm:344-352, 404-441`). Every allele of an input
//! record carries the record's whole list; a renderer keeps the records whose
//! `allele` is the row's.

use std::sync::Arc;

use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// One source's annotation of one variant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct CustomAnnotations {
    /// The source's `short_name`, the output key.
    pub short_name: Arc<str>,
    /// The matched records in file order: at most `num_records`, then one named
    /// `...`.
    pub records: Vec<CustomRecord>,
    /// Whether the source matched a record (`$record_count`, `File.pm:339`):
    /// a source asked for statistics leaves an entry when it matched none,
    /// and one whose every match lent its fields to an earlier source's
    /// namesake (`:405-415`) lists no record of its own.
    pub matched: bool,
    /// The summary statistics; `None` when the source computes none.
    pub stats: Option<CustomStats>,
}

impl CustomAnnotations {
    /// One source's annotation of one variant, from each of its fields.
    pub fn new(
        short_name: Arc<str>,
        records: Vec<CustomRecord>,
        matched: bool,
        stats: Option<CustomStats>,
    ) -> Self {
        Self {
            short_name,
            records,
            matched,
            stats,
        }
    }
}

/// The summary statistics of one source over one variant
/// (`File.pm:344-352`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum CustomStats {
    /// No statistic was computed although statistics were asked for: `sum`
    /// alone of `sum`, `min` and `max` was named and no matched record
    /// carried a score, so Ensembl VEP's hash is never created.
    Undefined,
    /// The computed statistics by name (`min`, `max`, `mean`, `sum`, `count`),
    /// in the order the source computed them; empty when none was computed.
    Values(IndexMap<String, CustomValue>),
}

/// One matched record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct CustomRecord {
    /// The record's identifier, its coordinates under `coords`, or `...`.
    pub name: String,
    /// The input allele an `exact` VCF record matched; `None` for a record that
    /// applies to every allele of the input record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allele: Option<Vec<u8>>,
    /// The values of the `fields` the source names, by field.
    #[serde(default, skip_serializing_if = "IndexMap::is_empty")]
    pub fields: IndexMap<String, CustomFieldValue>,
    /// The record's score, which the statistics run over; `None` when none was
    /// read.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<CustomValue>,
}

impl CustomRecord {
    /// One matched record, from each of its fields.
    pub fn new(
        name: String,
        allele: Option<Vec<u8>>,
        fields: IndexMap<String, CustomFieldValue>,
        score: Option<CustomValue>,
    ) -> Self {
        Self {
            name,
            allele,
            fields,
            score,
        }
    }
}

/// The value of one field of a record: the column's text, the list a
/// `FILTER` column splits into on `;` (`File/VCF.pm:222`), or a value the
/// record leaves undefined (a GFF window's `feature_id` when the record has no
/// `ID` attribute, `File/GFF.pm:199`), which the text formats print empty and
/// the JSON format `null`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum CustomFieldValue {
    Text(String),
    List(Vec<String>),
    Undefined,
}

/// A score or statistic as the source read or computed it: the text of a
/// file's column, a number, or a value the file leaves undefined (a `.`
/// QUAL).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum CustomValue {
    Text(String),
    Number(f64),
    Undefined,
}
