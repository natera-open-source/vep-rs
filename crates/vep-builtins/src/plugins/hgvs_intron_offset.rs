// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVSIntronOffset plugin: the intron offsets of the HGVSc coordinates, as
//! Ensembl VEP's `HGVSIntronOffset.pm` writes them.
//!
//! **Usage:** `--hgvs --plugin HGVSIntronOffset` (the plugin takes no parameter;
//! a positional one is not read, and is logged; a `key=value` one is refused; a
//! run without `--hgvs` is refused at load)
//!
//! Output fields: `HGVS_IntronStartOffset` and `HGVS_IntronEndOffset`, the
//! signed distance after the exon coordinate of the HGVSc start and end
//! positions as `hgvs_transcript` keeps them (`:64-65`;
//! TranscriptVariationAllele.pm:1530-1545): `45` for `c.123+45`, `-10` for
//! `c.124-10`, `0` for an exonic position, read before the two positions are
//! ordered for printing. A row without an HGVSc carries nothing (`:62`), a
//! structural variant's rows included (`BaseVepPlugin.pm:100`).
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`]
//! via the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use indexmap::IndexMap;

use crate::params::parse_none;
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// The start offset field (`HGVSIntronOffset.pm:54`).
const START_FIELD: &str = "HGVS_IntronStartOffset";

/// The end offset field (`HGVSIntronOffset.pm:55`).
const END_FIELD: &str = "HGVS_IntronEndOffset";

/// HGVS intron offset plugin.
#[derive(Debug, Default)]
pub struct HgvsIntronOffsetPlugin;

impl HgvsIntronOffsetPlugin {
    /// Create the plugin; `--plugin HGVSIntronOffset` carries no parameter to
    /// initialize.
    pub fn new() -> Self {
        Self
    }
}

impl BuiltinPlugin for HgvsIntronOffsetPlugin {
    fn name(&self) -> &str {
        "HGVSIntronOffset"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                START_FIELD.to_string(),
                "HGVS intron start offset".to_string(),
            ),
            (END_FIELD.to_string(), "HGVS intron end offset".to_string()),
        ]
    }

    /// A positional parameter is not read, as the Ensembl VEP module reads
    /// none (HGVSIntronOffset.pm defines no `new`); it is logged. A `key=value`
    /// list is refused as every built-in refuses an undeclared key.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        parse_none(self.name(), params)
    }

    fn requires_hgvs(&self) -> bool {
        true
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if variant.is_structural || consequence.hgvsc.is_none() {
            return Ok(result);
        }
        if let Some((start, end)) = consequence.hgvs_intron_offsets {
            result.insert(START_FIELD.to_string(), start.to_string());
            result.insert(END_FIELD.to_string(), end.to_string());
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::FeatureType;

    fn row(hgvsc: Option<&str>, offsets: Option<(i32, i32)>) -> TranscriptConsequence {
        TranscriptConsequence {
            hgvsc: hgvsc.map(str::to_string),
            hgvs_intron_offsets: offsets,
            feature_type: FeatureType::Transcript,
            ..TranscriptConsequence::default()
        }
    }

    fn snv() -> InputVariant {
        InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"T".to_vec())
    }

    fn offsets(
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Vec<(String, String)> {
        HgvsIntronOffsetPlugin::new()
            .run(consequence, variant)
            .unwrap()
            .into_iter()
            .collect()
    }

    fn pair(start: &str, end: &str) -> Vec<(String, String)> {
        vec![
            (START_FIELD.to_string(), start.to_string()),
            (END_FIELD.to_string(), end.to_string()),
        ]
    }

    /// Ensembl VEP writes `HGVS_IntronStartOffset=45;HGVS_IntronEndOffset=-10`
    /// for `c.123+45_124-10del`, `0` and `0` for an exonic `c.300A>G`, and the
    /// negative distances of an acceptor-side position.
    #[test]
    fn the_two_offsets_are_written_signed_with_zero_for_exonic_positions() {
        let intronic = row(Some("ENST1.1:c.123+45_124-10del"), Some((45, -10)));
        assert_eq!(offsets(&intronic, &snv()), pair("45", "-10"));
        let exonic = row(Some("ENST1.1:c.300A>G"), Some((0, 0)));
        assert_eq!(offsets(&exonic, &snv()), pair("0", "0"));
        let acceptor = row(Some("ENST1.1:c.124-3C>T"), Some((-3, -3)));
        assert_eq!(offsets(&acceptor, &snv()), pair("-3", "-3"));
    }

    #[test]
    fn a_row_without_hgvsc_carries_nothing() {
        assert!(offsets(&row(None, None), &snv()).is_empty());
        assert!(offsets(&row(None, Some((0, 0))), &snv()).is_empty());
    }

    #[test]
    fn a_structural_variant_row_carries_nothing() {
        let mut deletion = InputVariant::new("21".into(), 100, 500, b"N".to_vec(), b"-".to_vec());
        deletion.is_structural = true;
        let tc = row(Some("ENST1.1:c.123+45_124-10del"), Some((45, -10)));
        assert!(offsets(&tc, &deletion).is_empty());
    }

    /// The plugin needs the HGVS notations; it reads no parameter, as the
    /// Ensembl VEP module reads none (HGVSIntronOffset.pm defines no `new`), so
    /// a positional one is not read and the plugin loads; a `key=value` list is
    /// refused as the undeclared key every built-in refuses.
    #[test]
    fn the_plugin_needs_the_hgvs_notations_and_reads_no_parameter() {
        let mut plugin = HgvsIntronOffsetPlugin::new();
        assert!(plugin.requires_hgvs());
        assert!(plugin.init(&[]).is_ok());
        assert!(plugin.init(&["x".to_string()]).is_ok());
        assert_eq!(
            plugin
                .init(&["file=x".to_string()])
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: HGVSIntronOffset: unknown parameter 'file'; \
             HGVSIntronOffset takes no parameter"
        );
    }
}
