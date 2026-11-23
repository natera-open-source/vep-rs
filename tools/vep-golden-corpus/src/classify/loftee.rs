// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The LOFTEE field class of a GERP walk without the conservation database: the plugin
//! Ensembl VEP runs on GRCh37 dies where its walk crosses a whole exon, so the row carries
//! none of the plugin's fields, where vep-rs totals the exon's bases and writes them.

use super::fields::{FieldModels, FieldQuestion, FieldRow, Verdict};

/// The class of `docs/divergences/loftee-gerp-walk-without-the-conservation-database.md`:
/// on the `master` branch of the LOFTEE repository gerp_dist.pl:59-60 read the GERP total of
/// each whole exon the walk crosses from the conservation database and :129 call `prepare`
/// on its handle, which LoF.pm:119-120 leave as the text `false` in a run without
/// `conservation_file`.
pub const LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE: &str =
    "loftee_gerp_walk_without_conservation_database";

/// The plugin's four fields: the ones the class takes, and the ones it reads of both rows
/// whichever of them the run compares.
pub const LOFTEE_FIELDS: [&str; 4] = ["LoF", "LoF_filter", "LoF_flags", "LoF_info"];

/// The terms of a row whose GERP walk the plugin starts (LoF.pm:180 and :251).
const TRUNCATING_TERMS: [&str; 2] = ["frameshift_variant", "stop_gained"];

/// The `LoF_info` item that holds the walk's GERP total.
const GERP_DIST_ITEM: &str = "GERP_DIST";

/// Whether a value is none: the empty text, or the `-` of a field a row lacks.
fn unwritten(value: &str) -> bool {
    matches!(value, "" | "-")
}

/// The value of `name` on `row`, `-` for a field the row was not read with.
fn field<'a>(row: &'a FieldRow, name: &str) -> &'a str {
    row.get(name).map_or("-", String::as_str)
}

/// Whether a `LoF_info` value holds the walk's GERP total among its `KEY:VALUE` items, read
/// at either separator a format joins them with (`,`, and `&` in the VCF format), as the
/// format comparator reads the value.
fn holds_the_gerp_total(info: &str) -> bool {
    info.split([',', '&'])
        .any(|item| item.split(':').next() == Some(GERP_DIST_ITEM))
}

/// `loftee_gerp_walk_without_conservation_database`: one of the plugin's four fields that
/// Ensembl VEP leaves unwritten and vep-rs writes, on a key whose consequence set holds
/// `stop_gained` or `frameshift_variant`, where Ensembl VEP's row carries none of the four
/// and vep-rs's row carries `LoF` beside a `LoF_info` that holds a `GERP_DIST` item, or
/// beside no `LoF_info`, as a run under `skip_lof_info` writes the row (LoF.pm:243 and
/// :339). A question without both rows passes.
pub fn gerp_walk_without_conservation_database(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let (Some(perl_row), Some(rust_row)) = (question.perl_row, question.rust_row) else {
        return Verdict::Pass;
    };
    let truncating = question
        .key
        .consequences
        .split(',')
        .any(|term| TRUNCATING_TERMS.contains(&term));
    let info = field(rust_row, "LoF_info");
    if LOFTEE_FIELDS.contains(&question.name)
        && truncating
        && unwritten(question.vep)
        && !unwritten(question.vep_rs)
        && LOFTEE_FIELDS
            .iter()
            .all(|name| unwritten(field(perl_row, name)))
        && !unwritten(field(rust_row, "LoF"))
        && (unwritten(info) || holds_the_gerp_total(info))
    {
        return Verdict::Class(LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "loftee_tests.rs"]
mod tests;
