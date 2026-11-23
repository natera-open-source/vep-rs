// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a `<CNV:TR>` allele read symbolically whose record states fewer bases than
//! the run it spans, or exactly as many: Ensembl VEP names `feature_elongation` for every
//! tandem repeat (Utils/VariationEffect.pm:246-256 and :337-348), and vep-rs names
//! `feature_truncation` for the first and neither term for the second. Every symbolic
//! `<CNV:TR>` allele over one run writes the Allele `tandem_repeat` at one Location, so
//! the rows of a gain, of a loss and of an allele as long as the run share a key, and an
//! entry names the records whose rows differ ([`records_of_both`]).

use super::keys::PositionRow;
use std::collections::BTreeSet;

/// The class of `docs/divergences/cnv-tr-symbolic-loss-named-a-gain.md`.
pub const CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN: &str = "cnv_tr_symbolic_loss_named_a_gain";

/// The Allele both engines write on the row of a `<CNV:TR>` allele read symbolically.
const TANDEM_REPEAT_ALLELE: &str = "tandem_repeat";

/// The term Ensembl VEP alone writes, and the term vep-rs alone writes in its place for an
/// allele that states a loss.
const REFERENCE_ONLY_TERM: &str = "feature_elongation";
const CANDIDATE_ONLY_TERM: &str = "feature_truncation";

/// The names, each an `Uploaded_variation`, of the records that write a row among
/// `reference_rows` and a row among `candidate_rows`: the records on which the set of the
/// one stands against the set of the other on a key several records share.
pub fn records_of_both<'a>(
    reference_rows: &'a [PositionRow],
    candidate_rows: &'a [PositionRow],
) -> BTreeSet<&'a String> {
    let reference: BTreeSet<&String> = reference_rows
        .iter()
        .map(|row| &row.uploaded_variation)
        .collect();
    let candidate = candidate_rows.iter().map(|row| &row.uploaded_variation);
    candidate.filter(|name| reference.contains(name)).collect()
}

/// `cnv_tr_symbolic_loss_named_a_gain` for a key whose Ensembl VEP set `perl_csq` and
/// vep-rs set `rust_csq` differ as the class documents, read from the key and the two sets
/// alone: the key is a Transcript's, its `allele` is `tandem_repeat`, Ensembl VEP's set less
/// vep-rs's is `feature_elongation`, and vep-rs's less Ensembl VEP's is `feature_truncation`
/// (an allele that states a loss) or nothing, vep-rs's set then holding no
/// `feature_truncation` (an allele as long as its run). A key that misses one of these takes
/// no class.
pub fn classify_symbolic_loss_pair(
    feature_type: &str,
    allele: &str,
    perl_csq: &str,
    rust_csq: &str,
) -> Option<&'static str> {
    if feature_type != "Transcript" || allele != TANDEM_REPEAT_ALLELE {
        return None;
    }
    let perl: BTreeSet<&str> = perl_csq.split(',').collect();
    let rust: BTreeSet<&str> = rust_csq.split(',').collect();
    let perl_only: Vec<&str> = perl.difference(&rust).copied().collect();
    let rust_only: Vec<&str> = rust.difference(&perl).copied().collect();
    let names_a_loss = rust_only == [CANDIDATE_ONLY_TERM];
    let names_neither = rust_only.is_empty() && !rust.contains(CANDIDATE_ONLY_TERM);
    (perl_only == [REFERENCE_ONLY_TERM] && (names_a_loss || names_neither))
        .then_some(CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN)
}

#[cfg(test)]
#[path = "tandem_repeat_loss_tests.rs"]
mod tests;
