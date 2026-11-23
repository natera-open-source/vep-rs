// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reports: the JSON, the markdown, the two tables of one-sided tuples and the release
//! block.
//!
//! The JSON report has the layout of CPython's `json.dump(report, fh, indent=2)` followed by
//! one newline: keys in the order they were set, ASCII-only strings, a ratio rounded to six
//! decimals and written as `float.__repr__` writes it.

use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::failure::Failure;
use crate::masks::{
    masks_by_release, unknown_release, MaskResult, NamedCounts, MASKS_RETIRED_IN_116_2,
};
use crate::metrics::{DiscordantRecord, Metrics};
use crate::pyfmt::{push_json_str, py_float_repr};
use crate::registry::{
    check_every_rule_has_a_mask, excludes_with_cache, one_sided_classes_116_2,
    BND_BARE_REFERENCE_ALLELE_ROW_BUCKET, BND_MATE_32_BIT_BUFFER_BUCKET, BND_OWN_ROWS_LOST_BUCKET,
    BND_SYNONYM_MATE_ROWS_LOST_BUCKET, CNV_TR_EMPTY_EXPANSION_BUCKET,
    CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET, DOT_BREAKEND_PAIR_BUCKET,
    FASTA_NAMED_SLICE_INTERGENIC_BUCKET, GVCF_VARIANT_SITE_DROPPED_BUCKET,
    MISSING_ALT_ANNOTATED_BUCKET, NON_INTEGER_END_BUCKET, NON_NUMERIC_POS_NUMIFIED_BUCKET,
    REFERENCE_SKIPPED_RECORD_BUCKET, STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
};
use crate::rules::table::DivergenceClass;
use crate::strings::{Strings, TupleSet};
use crate::text::{os_failure, py_join};

/// Everything the four report files are written from. Lists keyed by a file's name are in
/// the order the files were compared; `per_type` is in the order of the type names.
pub struct Reports<'a> {
    pub strings: &'a Strings,
    pub assembly: &'a str,
    pub reference_release: &'a str,
    pub scored_engine: &'a str,
    /// Whether a cache authority was given and read.
    pub cache_supplied: bool,
    pub file_metrics: &'a [(String, Metrics)],
    pub per_type: &'a [(String, Metrics)],
    pub overall: &'a Metrics,
    pub adjusted_overall: &'a Metrics,
    pub adjusted_file_metrics: &'a [(String, Metrics)],
    /// The masks over the union of all files' tuples.
    pub masks: &'a MaskResult,
    pub masks_by_file: &'a [(String, MaskResult)],
    /// The size of the tandem-repeat class before its mask, `(rust, perl)`, over the union.
    pub cnvtr_total: (u64, u64),
    pub cnvtr_total_by_file: &'a [(String, (u64, u64))],
    pub discordants: &'a [DiscordantRecord],
    pub adjusted_discordants: &'a [DiscordantRecord],
}

/// A JSON value. An integer is a count of up to 64 bits or a signed difference of two.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i128),
    Float(f64),
    Str(String),
    Array(Vec<Json>),
    Object(JsonObject),
}

impl From<bool> for Json {
    fn from(value: bool) -> Json {
        Json::Bool(value)
    }
}

impl From<u64> for Json {
    fn from(value: u64) -> Json {
        Json::Int(i128::from(value))
    }
}

impl From<i64> for Json {
    fn from(value: i64) -> Json {
        Json::Int(i128::from(value))
    }
}

impl From<&str> for Json {
    fn from(value: &str) -> Json {
        Json::Str(value.to_string())
    }
}

impl From<JsonObject> for Json {
    fn from(value: JsonObject) -> Json {
        Json::Object(value)
    }
}

impl<T: Into<Json>> From<Option<T>> for Json {
    fn from(value: Option<T>) -> Json {
        value.map_or(Json::Null, Into::into)
    }
}

impl Json {
    /// The value as CPython's `json.dumps(value, indent=2)` writes it: no trailing newline.
    pub fn dumps(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out
    }

    fn write(&self, out: &mut String, level: usize) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(value) => out.push_str(if *value { "true" } else { "false" }),
            Json::Int(value) => out.push_str(&value.to_string()),
            Json::Float(value) if value.is_nan() => out.push_str("NaN"),
            Json::Float(value) if value.is_infinite() => {
                out.push_str(if *value < 0.0 {
                    "-Infinity"
                } else {
                    "Infinity"
                });
            }
            Json::Float(value) => out.push_str(&py_float_repr(*value)),
            Json::Str(value) => push_json_str(out, value),
            Json::Array(items) if items.is_empty() => out.push_str("[]"),
            Json::Object(object) if object.entries.is_empty() => out.push_str("{}"),
            Json::Array(items) => {
                out.push('[');
                for (n, item) in items.iter().enumerate() {
                    out.push_str(if n == 0 { "\n" } else { ",\n" });
                    out.push_str(&"  ".repeat(level + 1));
                    item.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&"  ".repeat(level));
                out.push(']');
            }
            Json::Object(object) => {
                out.push('{');
                for (n, (key, value)) in object.entries.iter().enumerate() {
                    out.push_str(if n == 0 { "\n" } else { ",\n" });
                    out.push_str(&"  ".repeat(level + 1));
                    push_json_str(out, key);
                    out.push_str(": ");
                    value.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&"  ".repeat(level));
                out.push('}');
            }
        }
    }
}

/// A JSON object with the order and the replacement rule of a Python `dict`: keys in the
/// order they were first set, a key set again keeping its place and taking the new value.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct JsonObject {
    entries: Vec<(String, Json)>,
}

impl JsonObject {
    pub fn new() -> JsonObject {
        JsonObject::default()
    }

    pub fn set(&mut self, key: &str, value: impl Into<Json>) {
        let value = value.into();
        match self.entries.iter_mut().find(|(held, _)| held == key) {
            Some(entry) => entry.1 = value,
            None => self.entries.push((key.to_string(), value)),
        }
    }

    /// Sets every entry of `other`, in its order.
    pub fn update(&mut self, other: JsonObject) {
        for (key, value) in other.entries {
            self.set(&key, value);
        }
    }

    pub fn get(&self, key: &str) -> Option<&Json> {
        self.entries
            .iter()
            .find(|(held, _)| held == key)
            .map(|(_, value)| value)
    }

    /// The keys in order.
    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(key, _)| key.as_str())
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// `round(x, 6)`: the f64 nearest to the half-to-even six-decimal rounding of x's exact
/// binary value.
pub fn py_round6(x: f64) -> f64 {
    format!("{x:.6}").parse().expect("fixed-point float text")
}

/// `f"{n:,}"`.
pub fn py_grouped(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn count(set: &TupleSet) -> u64 {
    set.len() as u64
}

/// The entries of a list keyed by a file's name as a by-file map of the report holds them:
/// by name in bytewise order, a name listed twice standing once with its last value.
fn by_name<T>(entries: &[(String, T)]) -> BTreeMap<&str, &T> {
    entries
        .iter()
        .map(|(name, value)| (name.as_str(), value))
        .collect()
}

fn named_counts(counts: &NamedCounts) -> JsonObject {
    let mut object = JsonObject::new();
    for (name, n) in counts {
        object.set(name, *n);
    }
    object
}

/// One row of counts and ratios as the report states it, the ratios rounded to six
/// decimals; a pooled row adds its five union-based counts.
pub fn metrics_to_dict(m: &Metrics) -> JsonObject {
    let mut d = JsonObject::new();
    d.set("total_input", m.total_input);
    d.set("perl_variants", m.perl_variants);
    d.set("rust_variants", m.rust_variants);
    d.set("perl_tuples", m.perl_tuples);
    d.set("rust_tuples", m.rust_tuples);
    d.set("intersection", m.intersection);
    d.set("only_perl", m.only_perl);
    d.set("only_rust", m.only_rust);
    d.set("precision", Json::Float(py_round6(m.precision)));
    d.set("recall", Json::Float(py_round6(m.recall)));
    d.set("f1", Json::Float(py_round6(m.f1)));
    if m.perl_tuples_union.is_some() {
        d.set("perl_tuples_union", m.perl_tuples_union);
        d.set("rust_tuples_union", m.rust_tuples_union);
        d.set("intersection_union", m.intersection_union);
        d.set("cross_file_duplicate_perl", m.cross_file_duplicate_perl);
        d.set("cross_file_duplicate_rust", m.cross_file_duplicate_rust);
    }
    d
}

fn metrics_by_name(rows: &[(String, Metrics)]) -> JsonObject {
    let mut object = JsonObject::new();
    for (name, m) in rows {
        object.set(name, metrics_to_dict(m));
    }
    object
}

/// A by-file map of one count a file: the files whose count is not zero, by name.
fn counts_by_file(
    per_file: &BTreeMap<&str, &MaskResult>,
    of: impl Fn(&MaskResult) -> u64,
) -> JsonObject {
    let mut object = JsonObject::new();
    for (name, masks) in per_file {
        let n = of(masks);
        if n != 0 {
            object.set(name, n);
        }
    }
    object
}

/// A by-file map of `{"rust": .., "perl": ..}` followed by the named counts of the same
/// mask: the files whose `rust` or `perl` entry is not zero, by name.
fn sides_by_file<'a, T>(
    per_file: &BTreeMap<&str, &'a T>,
    of: impl Fn(&'a T) -> (u64, u64, Option<&'a NamedCounts>),
) -> JsonObject {
    let mut object = JsonObject::new();
    for (name, value) in per_file {
        let (rust, perl, named) = of(value);
        let mut entry = JsonObject::new();
        entry.set("rust", rust);
        entry.set("perl", perl);
        if let Some(named) = named {
            entry.update(named_counts(named));
        }
        let charged = |side: &str| !matches!(entry.get(side), Some(Json::Int(0)));
        if charged("rust") || charged("perl") {
            object.set(name, entry);
        }
    }
    object
}

/// What the report of a release after 115.2 adds to its `adjusted` block: the release, the
/// masks applied and the 115.2 masks retired with the reason, the exclusions of each 116.2
/// class over the union and per file, whether a cache authority was supplied, and every
/// one-sided class of the registry with its definition, this comparator's verdict and its
/// count. A registry rule this comparator has no mask for is refused.
pub fn release_report_block(
    overall: &MaskResult,
    per_file: &[(String, MaskResult)],
    cache_supplied: bool,
) -> Result<JsonObject, Failure> {
    release_block_of(one_sided_classes_116_2(), overall, per_file, cache_supplied)
}

fn release_block_of(
    rules: &[DivergenceClass],
    overall: &MaskResult,
    per_file: &[(String, MaskResult)],
    cache_supplied: bool,
) -> Result<JsonObject, Failure> {
    let excluded_by_bucket: [(&str, u64); 14] = [
        (
            REFERENCE_SKIPPED_RECORD_BUCKET,
            count(&overall.skipped_rust),
        ),
        (BND_OWN_ROWS_LOST_BUCKET, count(&overall.bnd_own_rust)),
        (
            BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
            count(&overall.bnd_synonym_rust),
        ),
        (FASTA_NAMED_SLICE_INTERGENIC_BUCKET, 0),
        (
            BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
            overall.bare_reference_rows,
        ),
        (BND_MATE_32_BIT_BUFFER_BUCKET, 0),
        (GVCF_VARIANT_SITE_DROPPED_BUCKET, count(&overall.gvcf_rust)),
        (
            DOT_BREAKEND_PAIR_BUCKET,
            count(&overall.dot_pair_rust) + count(&overall.dot_pair_perl),
        ),
        (
            MISSING_ALT_ANNOTATED_BUCKET,
            count(&overall.missing_alt_perl),
        ),
        (
            NON_INTEGER_END_BUCKET,
            count(&overall.non_integer_end_perl) + count(&overall.non_integer_end_rust),
        ),
        (
            CNV_TR_EMPTY_EXPANSION_BUCKET,
            count(&overall.cnvtr_empty_rust),
        ),
        (
            STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
            count(&overall.star_alone_perl),
        ),
        (
            NON_NUMERIC_POS_NUMIFIED_BUCKET,
            count(&overall.numified_pos_perl),
        ),
        (
            CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET,
            count(&overall.rn_zero_rust),
        ),
    ];
    let masked: Vec<&str> = excluded_by_bucket
        .iter()
        .map(|(bucket, _)| *bucket)
        .collect();
    check_every_rule_has_a_mask(rules, &masked)?;
    let Some(masks_applied) = masks_by_release(&overall.release) else {
        return Err(unknown_release(&overall.release));
    };
    let per_file = by_name(per_file);

    let mut block = JsonObject::new();
    block.set("reference_release", overall.release.as_str());
    block.set(
        "masks_applied",
        Json::Array(masks_applied.iter().map(|mask| Json::from(*mask)).collect()),
    );
    let mut retired = JsonObject::new();
    for (mask, reason) in MASKS_RETIRED_IN_116_2 {
        retired.set(mask, reason);
    }
    block.set("masks_not_applied", retired);
    block.set(
        "excluded_reference_skipped_rust",
        count(&overall.skipped_rust),
    );
    block.set("reference_skipped_records", overall.skipped_records);
    let mut skipped_by_file = JsonObject::new();
    for (name, masks) in &per_file {
        if !masks.skipped_rust.is_empty() {
            let mut entry = JsonObject::new();
            entry.set("rust", count(&masks.skipped_rust));
            entry.set("records", masks.skipped_records);
            skipped_by_file.set(name, entry);
        }
    }
    block.set("excluded_reference_skipped_by_file", skipped_by_file);
    block.set(
        "excluded_bnd_own_rows_lost_rust",
        count(&overall.bnd_own_rust),
    );
    block.set(
        "excluded_bnd_own_rows_lost_by_file",
        counts_by_file(&per_file, |masks| count(&masks.bnd_own_rust)),
    );
    block.set(
        "excluded_bnd_synonym_mate_rows_rust",
        count(&overall.bnd_synonym_rust),
    );
    block.set(
        "excluded_bnd_synonym_mate_rows_by_file",
        counts_by_file(&per_file, |masks| count(&masks.bnd_synonym_rust)),
    );
    block.set(
        "fasta_named_slice_intergenic_records",
        overall.fasta_intergenic_records,
    );
    block.set(
        "bnd_bare_reference_allele_rows",
        overall.bare_reference_rows,
    );
    block.set(
        "bnd_bare_reference_allele_rows_by_file",
        counts_by_file(&per_file, |masks| masks.bare_reference_rows),
    );
    block.set(
        "bnd_mate_32bit_buffer_records",
        overall.mate_overflow_records,
    );
    block.set(
        "bnd_mate_32bit_buffer_reference_tuples",
        overall.mate_overflow_perl_tuples,
    );
    block.set("excluded_gvcf_variant_site_rust", count(&overall.gvcf_rust));
    block.set("gvcf_variant_site_records", overall.gvcf_records);
    block.set(
        "excluded_dot_breakend_pair_rust",
        count(&overall.dot_pair_rust),
    );
    block.set(
        "excluded_dot_breakend_pair_perl",
        count(&overall.dot_pair_perl),
    );
    block.set(
        "excluded_missing_alt_annotated_perl",
        count(&overall.missing_alt_perl),
    );
    block.set(
        "excluded_non_integer_end_perl",
        count(&overall.non_integer_end_perl),
    );
    block.set(
        "excluded_non_integer_end_rust",
        count(&overall.non_integer_end_rust),
    );
    block.set(
        "excluded_cnv_tr_empty_expansion_rust",
        count(&overall.cnvtr_empty_rust),
    );
    block.set(
        "excluded_star_alone_svtype_perl",
        count(&overall.star_alone_perl),
    );
    block.set(
        "excluded_non_numeric_pos_perl",
        count(&overall.numified_pos_perl),
    );
    block.set("excluded_cnv_tr_rn_zero_rust", count(&overall.rn_zero_rust));
    block.set(
        "scored_engine_skip_route_residue_rows",
        count(&overall.scored_engine_skip_residue),
    );
    block.set(
        "scored_engine_skip_route_residue_by_file",
        counts_by_file(&per_file, |masks| count(&masks.scored_engine_skip_residue)),
    );
    block.set("breakend_masks_check_ran", cache_supplied);
    let mut classes = JsonObject::new();
    for rule in rules {
        let excluded = excluded_by_bucket
            .iter()
            .find(|(bucket, _)| *bucket == rule.bucket)
            .map_or(0, |(_, n)| *n);
        let mut class = JsonObject::new();
        class.set("definition", rule.definition);
        class.set("excludes", excludes_with_cache(rule.bucket));
        class.set("excluded", excluded);
        class.set("taxonomy_class", rule.taxonomy_class);
        class.set("perl_citation", rule.perl_citation);
        classes.set(rule.bucket, class);
    }
    block.set("divergence_classes", classes);
    Ok(block)
}

/// The sizes of the two unions the adjusted counts subtract: every scored-side layer and
/// every reference-side layer of the masks over the union of all files.
fn excluded_totals(masks: &MaskResult) -> (u64, u64) {
    (count(&masks.all_rust()), count(&masks.all_perl()))
}

/// The JSON report as its file holds it, the closing newline included. Its `adjusted` block
/// states, beside the adjusted row, every count the mask accounting is rebuilt from: the two
/// unions, each layer's own size, the overlaps of the transcript-selection layer, the class
/// total of the tandem-repeat class and the per-file maps; a release after 115.2 adds its
/// release block.
pub fn json_report(reports: &Reports<'_>) -> Result<String, Failure> {
    json_report_with(reports, excluded_totals(reports.masks))
}

fn json_report_with(reports: &Reports<'_>, excluded: (u64, u64)) -> Result<String, Failure> {
    let masks = reports.masks;
    let release_block = if reports.reference_release != "115.2" {
        Some(release_report_block(
            masks,
            reports.masks_by_file,
            reports.cache_supplied,
        )?)
    } else {
        None
    };
    let per_file = by_name(reports.masks_by_file);
    let (excluded_rust, excluded_perl) = excluded;

    let mut adjusted = metrics_to_dict(reports.adjusted_overall);
    adjusted.set("excluded_divergences", excluded_rust);
    adjusted.set("excluded_rust_extra", excluded_rust);
    adjusted.set("excluded_perl_extra", excluded_perl);
    adjusted.set("excluded_cnvtr_rust", count(&masks.cnvtr_rust));
    adjusted.set("excluded_cnvtr_perl", count(&masks.cnvtr_perl));
    adjusted.set(
        "excluded_transcript_selection_rust",
        count(&masks.tsel_rust),
    );
    adjusted.set(
        "excluded_transcript_selection_perl",
        count(&masks.tsel_perl),
    );
    adjusted.set(
        "transcript_selection_cnvtr_overlap_rust",
        masks.tsel_rust.intersection_len(&masks.cnvtr_rust) as u64,
    );
    adjusted.set("cnvtr_swap_pairs_total_rust", reports.cnvtr_total.0);
    adjusted.set("cnvtr_swap_pairs_total_perl", reports.cnvtr_total.1);
    adjusted.set("excluded_cross_chromosome_perl", count(&masks.xchr_perl));
    adjusted.set(
        "excluded_cross_chromosome_rust",
        count(&masks.xchr_orphan_rust),
    );
    adjusted.set(
        "excluded_transcript_selection_mate_rust",
        count(&masks.mate_rust),
    );
    adjusted.set("excluded_registry_swap_rust", count(&masks.swap_rust));
    adjusted.set("excluded_registry_swap_perl", count(&masks.swap_perl));
    adjusted.set(
        "excluded_registry_swap_by_bucket",
        named_counts(&masks.swap_by_bucket),
    );
    adjusted.set("excluded_non_ref_batch_rust", count(&masks.non_ref_rust));
    adjusted.set("excluded_non_ref_batch_perl", count(&masks.non_ref_perl));
    adjusted.set(
        "excluded_breakend_mate_local_read_rust",
        count(&masks.mate_read_rust),
    );
    adjusted.set(
        "excluded_breakend_mate_local_read_perl",
        count(&masks.mate_read_perl),
    );
    adjusted.set(
        "excluded_breakend_mate_local_read_by_kind",
        named_counts(&masks.mate_read_by_kind),
    );
    adjusted.set(
        "transcript_selection_cross_chromosome_overlap_perl",
        masks.tsel_perl.intersection_len(&masks.xchr_perl) as u64,
    );
    adjusted.set("cross_chromosome_check_ran", reports.cache_supplied);
    adjusted.set("scored_engine", reports.scored_engine);
    adjusted.set(
        "scored_engine_transcripts_absent_from_cache",
        masks.xchr_unresolvable_scored,
    );
    // Five of the per-file maps are left out of a report over no file; the other three are
    // written, empty.
    let any_file = !per_file.is_empty();
    if any_file {
        adjusted.set(
            "excluded_cross_chromosome_by_file",
            counts_by_file(&per_file, |masks| count(&masks.xchr_perl)),
        );
    }
    adjusted.set(
        "excluded_cnvtr_by_file",
        counts_by_file(&per_file, |masks| count(&masks.cnvtr_rust)),
    );
    adjusted.set(
        "excluded_transcript_selection_by_file",
        sides_by_file(&per_file, |masks: &MaskResult| {
            (count(&masks.tsel_rust), count(&masks.tsel_perl), None)
        }),
    );
    adjusted.set(
        "cnvtr_swap_pairs_total_by_file",
        sides_by_file(
            &by_name(reports.cnvtr_total_by_file),
            |total: &(u64, u64)| (total.0, total.1, None),
        ),
    );
    if any_file {
        adjusted.set(
            "excluded_cross_chromosome_rust_by_file",
            counts_by_file(&per_file, |masks| count(&masks.xchr_orphan_rust)),
        );
        adjusted.set(
            "excluded_transcript_selection_mate_by_file",
            counts_by_file(&per_file, |masks| count(&masks.mate_rust)),
        );
        adjusted.set(
            "excluded_registry_swap_by_file",
            sides_by_file(&per_file, |masks: &MaskResult| {
                (
                    count(&masks.swap_rust),
                    count(&masks.swap_perl),
                    Some(&masks.swap_by_bucket),
                )
            }),
        );
        adjusted.set(
            "excluded_non_ref_batch_by_file",
            sides_by_file(&per_file, |masks: &MaskResult| {
                (count(&masks.non_ref_rust), count(&masks.non_ref_perl), None)
            }),
        );
        adjusted.set(
            "excluded_breakend_mate_local_read_by_file",
            sides_by_file(&per_file, |masks: &MaskResult| {
                (
                    count(&masks.mate_read_rust),
                    count(&masks.mate_read_perl),
                    Some(&masks.mate_read_by_kind),
                )
            }),
        );
    }
    if let Some(release_block) = release_block {
        adjusted.update(release_block);
    }

    let mut report = JsonObject::new();
    report.set("overall", metrics_to_dict(reports.overall));
    report.set("per_file", metrics_by_name(reports.file_metrics));
    report.set("per_variant_type", metrics_by_name(reports.per_type));
    report.set("adjusted", adjusted);
    report.set(
        "per_file_adjusted",
        metrics_by_name(reports.adjusted_file_metrics),
    );
    let mut text = Json::Object(report).dumps();
    text.push('\n');
    Ok(text)
}

/// The line of a report after 115.2 that states each class of the release beside the
/// headline's total.
fn release_classes_line(release: &str, masks: &MaskResult) -> String {
    let n = |set: &TupleSet| py_grouped(count(set));
    format!(
        "- **{release} classes**: `{REFERENCE_SKIPPED_RECORD_BUCKET}` \
         {} vep-rs rows on {} records set aside; \
         `{BND_OWN_ROWS_LOST_BUCKET}` {} set aside; \
         `{BND_SYNONYM_MATE_ROWS_LOST_BUCKET}` {} set aside; \
         `{FASTA_NAMED_SLICE_INTERGENIC_BUCKET}` {} records counted, none set aside; \
         `{BND_BARE_REFERENCE_ALLELE_ROW_BUCKET}` {} reference-allele rows set aside; \
         `{BND_MATE_32_BIT_BUFFER_BUCKET}` {} records counted, none set aside; \
         `{GVCF_VARIANT_SITE_DROPPED_BUCKET}` {} vep-rs rows on {} sites set aside; \
         `{DOT_BREAKEND_PAIR_BUCKET}` {} reference and {} vep-rs rows set aside; \
         `{MISSING_ALT_ANNOTATED_BUCKET}` {} reference rows set aside; \
         `{NON_INTEGER_END_BUCKET}` {} reference and {} vep-rs rows set aside; \
         `{CNV_TR_EMPTY_EXPANSION_BUCKET}` {} vep-rs rows set aside; \
         `{STAR_ALONE_SVTYPE_ANNOTATED_BUCKET}` {} reference rows set aside; \
         `{NON_NUMERIC_POS_NUMIFIED_BUCKET}` {} reference rows set aside; \
         `{CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET}` {} vep-rs rows set aside; \
         `scored_engine_skip_route_residue` {} scored-engine rows \
         at Locations the reference wrote nothing for counted, none set aside",
        n(&masks.skipped_rust),
        py_grouped(masks.skipped_records),
        n(&masks.bnd_own_rust),
        n(&masks.bnd_synonym_rust),
        py_grouped(masks.fasta_intergenic_records),
        py_grouped(masks.bare_reference_rows),
        py_grouped(masks.mate_overflow_records),
        n(&masks.gvcf_rust),
        py_grouped(masks.gvcf_records),
        n(&masks.dot_pair_perl),
        n(&masks.dot_pair_rust),
        n(&masks.missing_alt_perl),
        n(&masks.non_integer_end_perl),
        n(&masks.non_integer_end_rust),
        n(&masks.cnvtr_empty_rust),
        n(&masks.star_alone_perl),
        n(&masks.numified_pos_perl),
        n(&masks.rn_zero_rust),
        n(&masks.scored_engine_skip_residue),
    )
}

/// The markdown report as its file holds it: the headline counts, the per-file table with
/// an adjusted F1 where it differs from the raw one, and the per-type table with a mark on
/// a type whose F1 is above zero and below 0.99. A release after 115.2 names itself, the
/// masks it applied and its classes.
pub fn markdown_report(reports: &Reports<'_>) -> String {
    markdown_report_with(reports, excluded_totals(reports.masks))
}

fn markdown_report_with(reports: &Reports<'_>, excluded: (u64, u64)) -> String {
    let masks = reports.masks;
    let overall = reports.overall;
    let (excluded_rust, excluded_perl) = excluded;
    let mut lines: Vec<String> = vec![
        "# SV Concordance Report".to_string(),
        String::new(),
        "## Overall".to_string(),
        String::new(),
    ];
    if reports.reference_release != "115.2" {
        let applied: Vec<String> = masks_by_release(reports.reference_release)
            .unwrap_or(&[])
            .iter()
            .map(|mask| format!("`{mask}`"))
            .collect();
        lines.push(format!(
            "- **Reference**: Ensembl VEP {}; masks applied: {}",
            reports.reference_release,
            applied.join(", ")
        ));
        lines.push(release_classes_line(reports.reference_release, masks));
    }
    if excluded_rust + excluded_perl > 0 {
        lines.push(format!(
            "- **F1 (adjusted)**: {:.4}  *(excluded {} Rust-extra + {} Perl-extra transcript divergences)*",
            reports.adjusted_overall.f1,
            py_grouped(excluded_rust),
            py_grouped(excluded_perl)
        ));
    }
    lines.push(format!("- **F1 (raw)**: {:.4}", overall.f1));
    lines.push(format!(
        "- **Total input variants**: {}",
        overall.total_input
    ));
    lines.push(format!("- **Perl tuples**: {}", overall.perl_tuples));
    lines.push(format!("- **Rust tuples**: {}", overall.rust_tuples));
    lines.push(format!("- **Intersection**: {}", overall.intersection));
    lines.push(format!("- **Only Perl**: {}", overall.only_perl));
    lines.push(format!("- **Only Rust**: {}", overall.only_rust));
    let per_file = by_name(reports.masks_by_file);
    let cnvtr_rust = count(&masks.cnvtr_rust);
    if cnvtr_rust > 0 {
        let files: Vec<&str> = per_file
            .iter()
            .filter(|(_, masks)| !masks.cnvtr_rust.is_empty())
            .map(|(name, _)| *name)
            .collect();
        let assembly = if reports.assembly.is_empty() {
            String::new()
        } else {
            format!(" {}", reports.assembly.to_uppercase())
        };
        // "N of M": the mask takes a part of the class, and N alone reads as the whole.
        let of_total = match reports.cnvtr_total.0 {
            0 => String::new(),
            total => format!(" of {}", py_grouped(total)),
        };
        lines.push(String::new());
        lines.push(format!(
            "SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion \
             class ({}{of_total} pairs{assembly}, {} only). Raw F1 unaffected.",
            py_grouped(cnvtr_rust),
            files.join(", ")
        ));
    }
    lines.push(String::new());

    lines.push("## Per-File Metrics".to_string());
    lines.push(String::new());
    lines.push(
        "| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |".to_string(),
    );
    lines.push(
        "|------|-------|-------------|-------------|-----------|--------|--------|".to_string(),
    );
    let adjusted_by_file = by_name(reports.adjusted_file_metrics);
    for (name, m) in by_name(reports.file_metrics) {
        let adjusted = adjusted_by_file.get(name).copied().unwrap_or(m);
        let adjusted_f1 = if adjusted.f1 != m.f1 {
            format!("{:.4}", adjusted.f1)
        } else {
            String::new()
        };
        lines.push(format!(
            "| {name} | {} | {} | {} | {} | {:.4} | {adjusted_f1} |",
            m.total_input, m.perl_tuples, m.rust_tuples, m.intersection, m.f1
        ));
    }
    lines.push(String::new());

    lines.push("## Per-Variant-Type Metrics".to_string());
    lines.push(String::new());
    lines.push(
        "| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |"
            .to_string(),
    );
    lines.push(
        "|--------------|-------|-------------|-------------|-----------|-----------|------|"
            .to_string(),
    );
    for (variant_type, m) in by_name(reports.per_type) {
        let flag = if m.f1 < 0.99 && m.f1 > 0.0 { " **" } else { "" };
        lines.push(format!(
            "| {variant_type} | {} | {} | {} | {} | {} | {:.4}{flag} |",
            m.total_input, m.perl_tuples, m.rust_tuples, m.only_perl, m.only_rust, m.f1
        ));
    }
    lines.push(String::new());
    lines.join("\n") + "\n"
}

/// Makes the file `path` and hands it to `write`. A file that cannot be made is a failure
/// that names it; one that cannot be written to is the system's error alone.
fn write_file(path: &Path, write: impl FnOnce(File) -> std::io::Result<()>) -> Result<(), Failure> {
    let file = File::create(path).map_err(|error| os_failure(&error, Some(path)))?;
    write(file).map_err(|error| os_failure(&error, None))
}

fn write_text(path: &Path, text: &str) -> Result<(), Failure> {
    write_file(path, |mut file| file.write_all(text.as_bytes()))
}

/// Writes a table of one-sided tuples: a header, then one tab-separated row a record in the
/// order given, the tuple's five columns as their texts.
pub fn write_discordant_tsv(
    path: &Path,
    strings: &Strings,
    records: &[DiscordantRecord],
) -> Result<(), Failure> {
    write_file(path, |file| {
        let mut out = BufWriter::new(file);
        out.write_all(
            b"source_file\tvariant_type\tdirection\tlocation\tallele\tfeature\tfeature_type\tconsequence\n",
        )?;
        for record in records {
            let [location, allele, feature, feature_type, consequence] =
                strings.resolve(record.tuple);
            writeln!(
                out,
                "{}\t{}\t{}\t{location}\t{allele}\t{feature}\t{feature_type}\t{consequence}",
                record.source_file,
                record.variant_type,
                record.direction.as_str()
            )?;
        }
        out.flush()
    })
}

/// Writes `concordance_report.json`, `concordance_report.md`, `discordant.tsv` and
/// `discordant_adjusted.tsv` into `output_dir`, a directory in Python's spelling: a file
/// under the directory `.` is named without `./`, as the failure to make it names it.
pub fn write_reports(output_dir: &Path, reports: &Reports<'_>) -> Result<(), Failure> {
    let excluded = excluded_totals(reports.masks);
    let json = json_report_with(reports, excluded)?;
    write_text(
        &py_join(output_dir, OsStr::new("concordance_report.json")),
        &json,
    )?;
    write_text(
        &py_join(output_dir, OsStr::new("concordance_report.md")),
        &markdown_report_with(reports, excluded),
    )?;
    write_discordant_tsv(
        &py_join(output_dir, OsStr::new("discordant.tsv")),
        reports.strings,
        reports.discordants,
    )?;
    write_discordant_tsv(
        &py_join(output_dir, OsStr::new("discordant_adjusted.tsv")),
        reports.strings,
        reports.adjusted_discordants,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;

    fn at<'a>(object: &'a JsonObject, key: &str) -> &'a Json {
        object
            .get(key)
            .unwrap_or_else(|| panic!("no key {key} in the object"))
    }

    fn object<'a>(object: &'a JsonObject, key: &str) -> &'a JsonObject {
        match at(object, key) {
            Json::Object(inner) => inner,
            other => panic!("{key} is {other:?}"),
        }
    }

    #[test]
    fn a_value_is_written_as_json_dump_with_an_indent_of_two_writes_it() {
        let mut nested = JsonObject::new();
        nested.set(
            "array",
            Json::Array(vec![
                Json::Int(1),
                Json::from("two"),
                Json::Array(vec![
                    Json::Float(3.5),
                    Json::Null,
                    Json::Bool(true),
                    Json::Bool(false),
                ]),
            ]),
        );
        nested.set("text", "\u{e9}\"\\\n\t\u{2028} \u{1f600} ~\u{7f}");
        let floats = [
            0.0,
            1.0,
            0.5,
            1e-05,
            5e-05,
            1e-06,
            0.1,
            123456789.125,
            1e16,
            1e22,
            0.30000000000000004,
        ];
        let mut value = JsonObject::new();
        value.set("empty_object", JsonObject::new());
        value.set("empty_array", Json::Array(Vec::new()));
        value.set("nested", nested);
        value.set(
            "floats",
            Json::Array(floats.iter().map(|x| Json::Float(*x)).collect()),
        );
        value.set("big", u64::MAX);
        let expected = r#"{
  "empty_object": {},
  "empty_array": [],
  "nested": {
    "array": [
      1,
      "two",
      [
        3.5,
        null,
        true,
        false
      ]
    ],
    "text": "\u00e9\"\\\n\t\u2028 \ud83d\ude00 ~\u007f"
  },
  "floats": [
    0.0,
    1.0,
    0.5,
    1e-05,
    5e-05,
    1e-06,
    0.1,
    123456789.125,
    1e+16,
    1e+22,
    0.30000000000000004
  ],
  "big": 18446744073709551615
}"#;
        assert_eq!(Json::Object(value).dumps(), expected);
        assert_eq!(Json::Float(f64::NAN).dumps(), "NaN");
        assert_eq!(Json::Float(f64::NEG_INFINITY).dumps(), "-Infinity");
    }

    #[test]
    fn a_key_set_again_keeps_its_place() {
        let mut object = JsonObject::new();
        object.set("rust", 1u64);
        object.set("perl", 2u64);
        let mut later = JsonObject::new();
        later.set("rust", 3u64);
        later.set("other", -4i64);
        object.update(later);
        assert_eq!(object.keys().collect::<Vec<_>>(), ["rust", "perl", "other"]);
        assert_eq!(object.get("rust"), Some(&Json::Int(3)));
        assert_eq!(
            Json::Object(object).dumps(),
            "{\n  \"rust\": 3,\n  \"perl\": 2,\n  \"other\": -4\n}"
        );
    }

    #[test]
    fn a_ratio_is_rounded_to_six_decimals_half_to_even_on_its_binary_value() {
        let cases = [
            (0.0, "0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.6666666666666666, "0.666667"),
            (0.1234565, "0.123456"),
            (0.1234575, "0.123457"),
            (5e-07, "0.0"),
            (1.5e-06, "2e-06"),
            (5e-05, "5e-05"),
            (1e-06, "1e-06"),
            (0.9999995, "1.0"),
            (0.99999949, "0.999999"),
            (1e-07, "0.0"),
            (0.8421052631578948, "0.842105"),
            (123.4567895, "123.456789"),
        ];
        for (x, expected) in cases {
            assert_eq!(py_float_repr(py_round6(x)), expected, "round({x:?}, 6)");
        }
    }

    #[test]
    fn a_count_is_grouped_in_thousands() {
        let cases = [
            (0, "0"),
            (7, "7"),
            (999, "999"),
            (1000, "1,000"),
            (12345, "12,345"),
            (1234567, "1,234,567"),
            (u64::MAX, "18,446,744,073,709,551,615"),
        ];
        for (n, expected) in cases {
            assert_eq!(py_grouped(n), expected);
        }
    }

    #[test]
    fn a_ratio_of_a_table_is_written_to_four_decimals() {
        let cases = [
            (0.99995, "1.0000"),
            (5e-05, "0.0001"),
            (0.12345, "0.1235"),
            (0.12355, "0.1235"),
            (0.6666666666666666, "0.6667"),
            (1.0, "1.0000"),
            (0.0, "0.0000"),
            (0.98999999, "0.9900"),
        ];
        for (x, expected) in cases {
            assert_eq!(format!("{x:.4}"), expected, "{x:?}");
        }
    }

    #[test]
    fn a_row_states_its_counts_and_rounded_ratios_and_a_pooled_row_its_union_counts() {
        let row = Metrics {
            total_input: 10,
            perl_variants: 8,
            rust_variants: 9,
            perl_tuples: 100,
            rust_tuples: 90,
            intersection: 80,
            only_perl: 20,
            only_rust: 10,
            precision: 0.8888888888888888,
            recall: 0.8,
            f1: 0.8421052631578948,
            ..Metrics::default()
        };
        assert_eq!(
            Json::Object(metrics_to_dict(&row)).dumps(),
            r#"{
  "total_input": 10,
  "perl_variants": 8,
  "rust_variants": 9,
  "perl_tuples": 100,
  "rust_tuples": 90,
  "intersection": 80,
  "only_perl": 20,
  "only_rust": 10,
  "precision": 0.888889,
  "recall": 0.8,
  "f1": 0.842105
}"#
        );
        let pooled = Metrics {
            precision: 0.0000005,
            recall: 0.00005,
            f1: 1.0,
            perl_tuples_union: Some(95),
            intersection_union: Some(77),
            cross_file_duplicate_perl: Some(5),
            ..row.clone()
        };
        assert_eq!(
            Json::Object(metrics_to_dict(&pooled)).dumps(),
            r#"{
  "total_input": 10,
  "perl_variants": 8,
  "rust_variants": 9,
  "perl_tuples": 100,
  "rust_tuples": 90,
  "intersection": 80,
  "only_perl": 20,
  "only_rust": 10,
  "precision": 0.0,
  "recall": 5e-05,
  "f1": 1.0,
  "perl_tuples_union": 95,
  "rust_tuples_union": null,
  "intersection_union": 77,
  "cross_file_duplicate_perl": 5,
  "cross_file_duplicate_rust": null
}"#
        );
        let fewer_in_the_files = Metrics {
            cross_file_duplicate_perl: Some(-3),
            ..pooled.clone()
        };
        assert_eq!(
            metrics_to_dict(&fewer_in_the_files).get("cross_file_duplicate_perl"),
            Some(&Json::Int(-3)),
            "a file sum below the union is written with its sign"
        );
        let without_the_first = Metrics {
            perl_tuples_union: None,
            ..pooled
        };
        assert_eq!(
            metrics_to_dict(&without_the_first).keys().last(),
            Some("f1"),
            "the five union counts stand or fall with the first"
        );
    }

    fn release_116_2() -> MaskResult {
        MaskResult {
            release: "116.2".to_string(),
            ..MaskResult::default()
        }
    }

    #[test]
    fn test_the_release_block_describes_every_116_2_rule_with_this_comparators_verdict() {
        let block = release_report_block(&release_116_2(), &[], true).unwrap();
        let classes = object(&block, "divergence_classes");
        let rules = one_sided_classes_116_2();
        assert_eq!(
            classes.keys().collect::<Vec<_>>(),
            rules.iter().map(|rule| rule.bucket).collect::<Vec<_>>()
        );
        for rule in rules {
            let class = object(classes, rule.bucket);
            assert_eq!(
                at(class, "excludes"),
                &Json::Bool(excludes_with_cache(rule.bucket)),
                "{}",
                rule.bucket
            );
            assert!(class.get("excluded").is_some(), "{}", rule.bucket);
        }
        let verdicts: Vec<&Json> = rules
            .iter()
            .map(|rule| at(object(classes, rule.bucket), "excludes"))
            .collect();
        let expected = [
            true, true, true, false, true, false, true, true, true, true, true, true, true, true,
        ]
        .map(Json::Bool);
        assert_eq!(verdicts, expected.iter().collect::<Vec<_>>());
        assert_eq!(at(&block, "bnd_bare_reference_allele_rows"), &Json::Int(0));
        assert_eq!(at(&block, "bnd_mate_32bit_buffer_records"), &Json::Int(0));
        assert_eq!(at(&block, "excluded_gvcf_variant_site_rust"), &Json::Int(0));
        let with_rows = MaskResult {
            bare_reference_rows: 3,
            ..release_116_2()
        };
        let with_rows = release_report_block(&with_rows, &[], true).unwrap();
        assert_eq!(
            at(&with_rows, "bnd_bare_reference_allele_rows"),
            &Json::Int(3)
        );
        assert_eq!(
            at(
                object(
                    object(&with_rows, "divergence_classes"),
                    "bnd_bare_reference_allele_row"
                ),
                "excluded"
            ),
            &Json::Int(3)
        );
        let mut with_a_sixth: Vec<DivergenceClass> = rules
            .iter()
            .map(|rule| DivergenceClass { ..*rule })
            .collect();
        with_a_sixth.push(DivergenceClass {
            bucket: "a_sixth_shape",
            definition: "A shape.",
            excludes: false,
            taxonomy_class: "x",
            perl_citation: "X.pm:1",
        });
        let refused = release_block_of(&with_a_sixth, &release_116_2(), &[], true).unwrap_err();
        assert_eq!(refused.kind, FailureKind::Assertion);
        assert!(refused.message.contains("a_sixth_shape"), "{refused}");
    }

    #[test]
    fn the_release_block_keeps_its_keys_in_order_and_states_the_cache() {
        let block = release_report_block(&release_116_2(), &[], false).unwrap();
        assert_eq!(
            block.keys().collect::<Vec<_>>(),
            [
                "reference_release",
                "masks_applied",
                "masks_not_applied",
                "excluded_reference_skipped_rust",
                "reference_skipped_records",
                "excluded_reference_skipped_by_file",
                "excluded_bnd_own_rows_lost_rust",
                "excluded_bnd_own_rows_lost_by_file",
                "excluded_bnd_synonym_mate_rows_rust",
                "excluded_bnd_synonym_mate_rows_by_file",
                "fasta_named_slice_intergenic_records",
                "bnd_bare_reference_allele_rows",
                "bnd_bare_reference_allele_rows_by_file",
                "bnd_mate_32bit_buffer_records",
                "bnd_mate_32bit_buffer_reference_tuples",
                "excluded_gvcf_variant_site_rust",
                "gvcf_variant_site_records",
                "excluded_dot_breakend_pair_rust",
                "excluded_dot_breakend_pair_perl",
                "excluded_missing_alt_annotated_perl",
                "excluded_non_integer_end_perl",
                "excluded_non_integer_end_rust",
                "excluded_cnv_tr_empty_expansion_rust",
                "excluded_star_alone_svtype_perl",
                "excluded_non_numeric_pos_perl",
                "excluded_cnv_tr_rn_zero_rust",
                "scored_engine_skip_route_residue_rows",
                "scored_engine_skip_route_residue_by_file",
                "breakend_masks_check_ran",
                "divergence_classes",
            ]
        );
        assert_eq!(at(&block, "breakend_masks_check_ran"), &Json::Bool(false));
        let unknown = MaskResult {
            release: "117.0".to_string(),
            ..MaskResult::default()
        };
        let refused = release_report_block(&unknown, &[], false).unwrap_err();
        assert_eq!(refused.kind, FailureKind::Value);
        assert_eq!(
            refused.message,
            "unknown reference release '117.0'; one of ('115.2', '116.2', '116.3')"
        );
        let unknown = MaskResult {
            release: "it's a\\b".to_string(),
            ..MaskResult::default()
        };
        let refused = release_report_block(&unknown, &[], false).unwrap_err();
        assert_eq!(
            refused.message,
            "unknown reference release \"it's a\\\\b\"; one of ('115.2', '116.2', '116.3')"
        );
    }
}
