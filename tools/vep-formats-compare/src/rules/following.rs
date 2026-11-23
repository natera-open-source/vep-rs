// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cells that follow a term pair.

use crate::consts::CONSEQUENCE_RANK;
use crate::consts::HGVSP_FIELDS;
use crate::entry::Entry;
use crate::fields::field_named;
use crate::pyfmt;
use crate::rules::term_pairs::IMPACT_FIELDS;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::FRAMESHIFT_TERM;
use crate::shared::STOP_RETAINED_TERM;
use std::collections::BTreeSet;

/// ensembl-variation Utils/Constants.pm %OVERLAP_CONSEQUENCES at release/116, SO term to impact.
pub const TERM_IMPACT: &[(&str, &str)] = &[
    ("transcript_ablation", "HIGH"),
    ("splice_acceptor_variant", "HIGH"),
    ("splice_donor_variant", "HIGH"),
    ("stop_gained", "HIGH"),
    ("frameshift_variant", "HIGH"),
    ("stop_lost", "HIGH"),
    ("start_lost", "HIGH"),
    ("transcript_amplification", "HIGH"),
    ("feature_elongation", "HIGH"),
    ("feature_truncation", "HIGH"),
    ("inframe_insertion", "MODERATE"),
    ("inframe_deletion", "MODERATE"),
    ("missense_variant", "MODERATE"),
    ("protein_altering_variant", "MODERATE"),
    ("TFBS_ablation", "MODERATE"),
    ("splice_donor_5th_base_variant", "LOW"),
    ("splice_region_variant", "LOW"),
    ("splice_donor_region_variant", "LOW"),
    ("splice_polypyrimidine_tract_variant", "LOW"),
    ("incomplete_terminal_codon_variant", "LOW"),
    ("start_retained_variant", "LOW"),
    ("stop_retained_variant", "LOW"),
    ("synonymous_variant", "LOW"),
    ("coding_sequence_variant", "MODIFIER"),
    ("mature_miRNA_variant", "MODIFIER"),
    ("5_prime_UTR_variant", "MODIFIER"),
    ("3_prime_UTR_variant", "MODIFIER"),
    ("non_coding_transcript_exon_variant", "MODIFIER"),
    ("intron_variant", "MODIFIER"),
    ("NMD_transcript_variant", "MODIFIER"),
    ("non_coding_transcript_variant", "MODIFIER"),
    ("coding_transcript_variant", "MODIFIER"),
    ("upstream_gene_variant", "MODIFIER"),
    ("downstream_gene_variant", "MODIFIER"),
    ("TFBS_amplification", "MODIFIER"),
    ("TF_binding_site_variant", "MODIFIER"),
    ("regulatory_region_ablation", "MODIFIER"),
    ("regulatory_region_amplification", "MODIFIER"),
    ("regulatory_region_variant", "MODIFIER"),
    ("intergenic_variant", "MODIFIER"),
    ("sequence_variant", "MODIFIER"),
];

/// The Downstream plugin's fields, written only on a row with `frameshift_variant` and no
/// term holding `splice` (Downstream.pm:97-99).
pub const DOWNSTREAM_FIELDS: &[&str] = &["DownstreamProtein", "ProteinLengthChange"];

/// The fields of NMD and LOFTEE, written only on a row with one of `TRUNCATING_TERMS`
/// (NMD.pm:65 and :83, LoF.pm:166-167).
pub const TRUNCATING_TERM_FIELDS: &[&str] = &["NMD", "LoF", "LoF_info", "LoF_filter", "LoF_flags"];

/// Membership only.
pub const TRUNCATING_TERMS: &[&str] = &[
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
];

/// The fields of REVEL and AlphaMissense, written only on a row with `MISSENSE_TERM`
/// (REVEL.pm:162, AlphaMissense.pm:242).
pub const MISSENSE_TERM_FIELDS: &[&str] = &["REVEL", "am_", "alphamissense"];

pub const MISSENSE_TERM: &str = "missense_variant";

/// The fields of ClinPred and of dbNSFP under its own prefix, written only on a row with
/// one of `PROTEIN_CHANGE_TERMS` (ClinPred.pm:69 and :106; dbNSFP.pm:157 and :290-291,
/// where the plugin filters by consequence).
pub const PROTEIN_CHANGE_TERM_FIELDS: &[&str] = &["ClinPred", "dbNSFP_"];

/// Membership only.
pub const PROTEIN_CHANGE_TERMS: &[&str] =
    &["missense_variant", "stop_lost", "stop_gained", "start_lost"];

/// The IMPACT a row of these terms carries, its most severe term's (OutputFactory.pm:1337
/// and :1344 read the first consequence by rank), or `None` for a term without a rank.
pub fn terms_impact(terms: &[&str]) -> Option<String> {
    let mut most_severe: Option<(i64, &str)> = None;
    for term in terms {
        let (_, rank) = CONSEQUENCE_RANK.iter().find(|(ranked, _)| ranked == term)?;
        if most_severe.is_none_or(|(least, _)| *rank < least) {
            most_severe = Some((*rank, term));
        }
    }
    let (_, term) = most_severe?;
    let impact = TERM_IMPACT
        .iter()
        .find(|(ranked, _)| *ranked == term)
        .map(|(_, impact)| *impact);
    // A ranked term without an impact ends the run as the `KeyError` of `TERM_IMPACT[term]`.
    Some(impact.map_or_else(
        || crate::uncaught!("{}", pyfmt::repr_str(term)),
        str::to_string,
    ))
}

/// Whether a row of these terms allows the plugin that writes the field `name` to run, or
/// `None` for a field no term decides: a field of Downstream, of NMD or LOFTEE, of REVEL or
/// AlphaMissense, or of ClinPred or dbNSFP.
pub fn terms_write_field(name: &str, terms: &[&str]) -> Option<bool> {
    if field_named(name, DOWNSTREAM_FIELDS) {
        return Some(
            terms.contains(&FRAMESHIFT_TERM) && !terms.iter().any(|term| term.contains("splice")),
        );
    }
    if field_named(name, TRUNCATING_TERM_FIELDS) {
        return Some(terms.iter().any(|term| TRUNCATING_TERMS.contains(term)));
    }
    if field_named(name, MISSENSE_TERM_FIELDS) {
        return Some(terms.contains(&MISSENSE_TERM));
    }
    if field_named(name, PROTEIN_CHANGE_TERM_FIELDS) {
        return Some(terms.iter().any(|term| PROTEIN_CHANGE_TERMS.contains(term)));
    }
    None
}

/// Whether an HGVSp is a frameshift's notation: `fs`, or `ext` where the frameshift
/// begins at the stop codon.
pub fn hgvsp_frameshift_form(value: &str) -> bool {
    let change = pyfmt::partition(value, ":p.").2;
    change.contains("fs") || change.contains("ext")
}

/// Whether the differing values `ev` and `av` of the field `name` on the reference
/// entry `e` and the candidate entry `a` of one key are the ones the two entries'
/// Consequence sets give: the IMPACT of each set, a field of a plugin a term gates
/// (`terms_write_field`) present exactly where the set allows the plugin to run, or an
/// HGVSp that is a frameshift's notation exactly on the side without
/// `stop_retained_variant`.
pub fn term_pair_following_cell(e: &Entry, a: &Entry, name: &str, ev: &str, av: &str) -> bool {
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    // A Consequence without a term is the one empty term.
    let expected_terms: BTreeSet<&str> = expected_set.split(',').collect();
    let actual_terms: BTreeSet<&str> = actual_set.split(',').collect();
    if expected_terms == actual_terms || ev == av {
        return false;
    }
    let expected_terms: Vec<&str> = expected_terms.into_iter().collect();
    let actual_terms: Vec<&str> = actual_terms.into_iter().collect();
    if field_named(name, IMPACT_FIELDS) {
        return terms_impact(&expected_terms).as_deref() == Some(ev)
            && terms_impact(&actual_terms).as_deref() == Some(av);
    }
    if let (Some(expected_written), Some(actual_written)) = (
        terms_write_field(name, &expected_terms),
        terms_write_field(name, &actual_terms),
    ) {
        return expected_written != actual_written
            && blank(ev) != expected_written
            && blank(av) != actual_written;
    }
    if field_named(name, HGVSP_FIELDS) {
        let expected_kept = expected_terms.contains(&STOP_RETAINED_TERM);
        let actual_kept = actual_terms.contains(&STOP_RETAINED_TERM);
        if expected_kept == actual_kept
            || blank(ev)
            || blank(av)
            || pyfmt::partition(ev, ":p.").0 != pyfmt::partition(av, ":p.").0
        {
            return false;
        }
        return hgvsp_frameshift_form(ev) != expected_kept
            && hgvsp_frameshift_form(av) != actual_kept;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plugin_field_is_written_where_the_terms_allow_its_plugin_to_run() {
        let shifted = ["frameshift_variant", "3_prime_UTR_variant"];
        let spliced = ["frameshift_variant", "splice_region_variant"];
        assert_eq!(terms_write_field("DownstreamProtein", &shifted), Some(true));
        assert_eq!(
            terms_write_field("proteinlengthchange", &spliced),
            Some(false)
        );
        assert_eq!(
            terms_write_field("NMD", &["splice_donor_variant"]),
            Some(true)
        );
        assert_eq!(
            terms_write_field("lof_info", &["stop_retained_variant"]),
            Some(false)
        );
        assert_eq!(terms_write_field("LOF_FLAGS", &["stop_gained"]), Some(true));
        assert_eq!(terms_write_field("LoF", &[]), Some(false));
        assert_eq!(terms_write_field("IMPACT", &shifted), None);
        assert_eq!(terms_write_field("LoF_", &shifted), None);
    }

    #[test]
    fn a_frameshift_notation_is_read_after_the_protein_name() {
        assert!(hgvsp_frameshift_form("ENSP1.1:p.Tyr22CysfsTer4"));
        assert!(hgvsp_frameshift_form("ENSP1.1:p.Ter23TyrextTer4"));
        assert!(hgvsp_frameshift_form("a:p.b:p.fs"));
        assert!(!hgvsp_frameshift_form("ENSP1.1:p.Ter23="));
        assert!(!hgvsp_frameshift_form("fs_ext.1:p.Tyr22_Ter23delinsCys"));
        assert!(!hgvsp_frameshift_form("Tyr22CysfsTer4"));
    }
}
