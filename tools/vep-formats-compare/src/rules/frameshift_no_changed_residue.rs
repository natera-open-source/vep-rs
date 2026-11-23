// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a frameshift description that names one residue on both sides of its
//! residue number.

use crate::pyfmt;
use crate::rules::hgvsp_stop::STOP_FRAMESHIFT_RE_PATTERN;
use crate::rules::whole_text;
use crate::shared::{blank, consequence_set, FRAMESHIFT_TERM};
use std::sync::LazyLock;

static FRAMESHIFT_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_FRAMESHIFT_RE_PATTERN));

/// The name a protein description gives a stop.
const STOP_NAME: &str = "Ter";

/// Whether a row's terms, joined by `,` or by `&`, hold `frameshift_variant`.
fn holds_frameshift(consequence: &str) -> bool {
    consequence_set(consequence)
        .split(',')
        .any(|held| held == FRAMESHIFT_TERM)
}

/// `hgvsp_frameshift_with_no_changed_residue` on a differing HGVSp: the candidate writes
/// none; the reference's description, a pair of prediction parentheses set aside, is a
/// frameshift with the unknown count whose two residue names are one name and no stop's;
/// and the reference row's terms (`consequence`) hold `frameshift_variant`, as do the
/// candidate row's (`candidate_consequence`) where the caller holds that row.
pub fn frameshift_with_no_changed_residue(
    consequence: &str,
    candidate_consequence: Option<&str>,
    ev: &str,
    av: &str,
) -> bool {
    if blank(ev) || !blank(av) {
        return false;
    }
    if !holds_frameshift(consequence) || !candidate_consequence.is_none_or(holds_frameshift) {
        return false;
    }
    let (_, separator, description) = pyfmt::partition(ev, ":p.");
    if separator.is_empty() {
        return false;
    }
    let description = description
        .strip_prefix('(')
        .and_then(|inner| inner.strip_suffix(')'))
        .unwrap_or(description);
    FRAMESHIFT_WHOLE
        .captures(description)
        .is_some_and(|named| &named[4] == "?" && named[1] == named[3] && &named[1] != STOP_NAME)
}
