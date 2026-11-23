// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! TSSDistance plugin: the distance from the transcription start site, as
//! Ensembl VEP's `TSSDistance.pm` writes it.
//!
//! **Usage:** `--plugin TSSDistance` or `--plugin TSSDistance,both_direction=1`
//!
//! Output field: `TSSDistance`. On a forward-strand transcript the distance is
//! the transcript's start minus the variant's end, on a reverse-strand one the
//! variant's start minus the transcript's end (`TSSDistance.pm:85-90`): positive
//! upstream of the transcription start site. An upstream row carries it
//! (`:93-97`); every other row carries it, zero or negative, under
//! `both_direction=1` alone (`:100-104`), the parameter read by its truth as
//! `params_to_hash` leaves it (`:56-59`), so `both_direction=0` is off. The
//! plugin writes on transcript rows of sequence and structural variants alike
//! (`variant_feature_types` names `BaseVariationFeature`, `:73-75`).
//!
//! A parameter without `=` leaves Ensembl VEP's hash empty and the flag off
//! (`BaseVepPlugin.pm:123`); it is not read here either, and is logged. A key
//! other than `both_direction` is refused.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`]
//! via the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use indexmap::IndexMap;

use crate::params::{warn_unread_parameters, PluginParams};
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "TSSDistance";

/// The one parameter (`TSSDistance.pm:59`).
const BOTH_DIRECTION: &str = "both_direction";

/// Transcription start site distance plugin.
#[derive(Debug, Default)]
pub struct TssDistancePlugin {
    both_direction: bool,
}

impl TssDistancePlugin {
    /// Create the plugin; `init` reads `both_direction`.
    pub fn new() -> Self {
        Self::default()
    }
}

impl crate::traits::sealed::Sealed for TssDistancePlugin {}

impl BuiltinPlugin for TssDistancePlugin {
    fn name(&self) -> &str {
        "TSSDistance"
    }

    /// The Ensembl VEP module declares `['BaseVariationFeature']`
    /// (TSSDistance.pm:73-75), the class a VariationFeature and a
    /// StructuralVariationFeature both inherit.
    fn annotates_structural_variants(&self) -> bool {
        true
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "Distance from the transcription start site".to_string(),
        )]
    }

    /// `both_direction=` by its truth (`:56-59`); a token without `=` is not
    /// read, as `params_to_hash` leaves it out of the hash (`BaseVepPlugin.pm:123`).
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let parsed = PluginParams::parse(self.name(), &[BOTH_DIRECTION], params)?;
        warn_unread_parameters(self.name(), parsed.positional(), "both_direction=1");
        self.both_direction = parsed.get(BOTH_DIRECTION).is_some_and(perl_true);
        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        let distance = tss_distance(consequence, variant);
        if distance > 0 || self.both_direction {
            result.insert(OUTPUT_FIELD.to_string(), distance.to_string());
        }
        Ok(result)
    }
}

/// `TSSDistance.pm:85-90`: the transcript's start minus the variant's end on
/// the forward strand, the variant's start minus the transcript's end on the
/// reverse strand. The row's feature span is the transcript's; the variant's
/// is its own, an allele of a `--minimal` split record measuring from its
/// trimmed span (the documented class `minimal_plugin_sees_the_rejoined_record`).
fn tss_distance(consequence: &TranscriptConsequence, variant: &InputVariant) -> i64 {
    let (transcript_start, transcript_end) = (
        consequence.feature_start as i64,
        consequence.feature_end as i64,
    );
    if consequence.strand == 1 {
        transcript_start - variant.end as i64
    } else {
        variant.start as i64 - transcript_end
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::FeatureType;

    /// A transcript row spanning `start..=end` on `strand`.
    fn row(strand: i8, start: u64, end: u64) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.strand = strand;
        consequence.feature_start = start;
        consequence.feature_end = end;
        consequence.feature_type = FeatureType::Transcript;
        consequence
    }

    fn span(start: u64, end: u64) -> InputVariant {
        InputVariant::new("21".into(), start, end, b"C".to_vec(), b"T".to_vec())
    }

    fn plugin(params: &[&str]) -> TssDistancePlugin {
        let mut plugin = TssDistancePlugin::new();
        plugin
            .init(&params.iter().map(|p| p.to_string()).collect::<Vec<_>>())
            .unwrap();
        plugin
    }

    fn distance(
        plugin: &TssDistancePlugin,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Option<String> {
        plugin
            .run(consequence, variant)
            .unwrap()
            .get(OUTPUT_FIELD)
            .cloned()
    }

    /// Ensembl VEP writes `TSSDistance=900` for a substitution at 9100 upstream
    /// of a forward transcript starting at 10000, and `TSSDistance=1500` for one
    /// at 21500 upstream of a reverse transcript ending at 20000.
    #[test]
    fn an_upstream_row_carries_the_distance_to_the_start_site() {
        let plugin = plugin(&[]);
        let forward = row(1, 10000, 20000);
        assert_eq!(
            distance(&plugin, &forward, &span(9100, 9100)).as_deref(),
            Some("900")
        );
        // The distance is measured from the variant's end on the forward strand.
        assert_eq!(
            distance(&plugin, &forward, &span(9000, 9100)).as_deref(),
            Some("900")
        );
        let reverse = row(-1, 10000, 20000);
        assert_eq!(
            distance(&plugin, &reverse, &span(21500, 21500)).as_deref(),
            Some("1500")
        );
        // ... and from the variant's start on the reverse strand.
        assert_eq!(
            distance(&plugin, &reverse, &span(21500, 21600)).as_deref(),
            Some("1500")
        );
    }

    /// An insertion between 9999 and 10000 is the pair `start 10000, end 9999`:
    /// one base upstream on the forward strand.
    #[test]
    fn an_insertion_measures_from_its_flanks() {
        let plugin = plugin(&[]);
        assert_eq!(
            distance(&plugin, &row(1, 10000, 20000), &span(10000, 9999)).as_deref(),
            Some("1")
        );
        assert_eq!(
            distance(&plugin, &row(-1, 10000, 20000), &span(20001, 20000)).as_deref(),
            Some("1")
        );
    }

    /// An allele of a `--minimal` split record measures from its own trimmed
    /// span, the one its row's alleles describe, as the same allele does when
    /// it is the record's only alternate: the `-` allele of `CAA/C/T` at
    /// 13826644 upstream of a reverse transcript ending at 13825635 measures
    /// 1010 from its own start 13826645, the insertion allele of `C/CA/T`
    /// there likewise, and on the forward strand the `T` allele of `CAA/C/TAA`
    /// at 9100 measures 900 from its own end. The record's span, the one its
    /// `Location` prints, would give 1009 and 898.
    #[test]
    fn a_split_allele_measures_from_its_own_span() {
        let plugin = plugin(&[]);
        let reverse = row(-1, 13824406, 13825635);
        let mut deletion = span(13826645, 13826646);
        deletion.original_start = Some(13826644);
        deletion.original_end = Some(13826646);
        assert_eq!(
            distance(&plugin, &reverse, &deletion).as_deref(),
            Some("1010")
        );
        let mut insertion = span(13826645, 13826644);
        insertion.original_start = Some(13826644);
        insertion.original_end = Some(13826644);
        assert_eq!(
            distance(&plugin, &reverse, &insertion).as_deref(),
            Some("1010")
        );
        assert_eq!(
            distance(&plugin, &reverse, &span(13826645, 13826646)).as_deref(),
            Some("1010")
        );
        let mut substitution = span(9100, 9100);
        substitution.original_start = Some(9100);
        substitution.original_end = Some(9102);
        assert_eq!(
            distance(&plugin, &row(1, 10000, 20000), &substitution).as_deref(),
            Some("900")
        );
    }

    #[test]
    fn rows_inside_or_downstream_carry_nothing_by_default() {
        let plugin = plugin(&[]);
        let forward = row(1, 10000, 20000);
        assert_eq!(distance(&plugin, &forward, &span(10000, 10000)), None);
        assert_eq!(distance(&plugin, &forward, &span(15000, 15000)), None);
        assert_eq!(distance(&plugin, &forward, &span(25000, 25000)), None);
        let reverse = row(-1, 10000, 20000);
        assert_eq!(distance(&plugin, &reverse, &span(20000, 20000)), None);
        assert_eq!(distance(&plugin, &reverse, &span(5000, 5000)), None);
    }

    /// Under `both_direction=1` Ensembl VEP writes the same arithmetic on every
    /// row: `0` at the start site, `-5000` inside, `-15000` downstream.
    #[test]
    fn both_direction_writes_zero_and_negative_distances() {
        let plugin = plugin(&["both_direction=1"]);
        let forward = row(1, 10000, 20000);
        assert_eq!(
            distance(&plugin, &forward, &span(10000, 10000)).as_deref(),
            Some("0")
        );
        assert_eq!(
            distance(&plugin, &forward, &span(15000, 15000)).as_deref(),
            Some("-5000")
        );
        assert_eq!(
            distance(&plugin, &forward, &span(25000, 25000)).as_deref(),
            Some("-15000")
        );
        assert_eq!(
            distance(&plugin, &row(-1, 10000, 20000), &span(5000, 5000)).as_deref(),
            Some("-15000")
        );
    }

    /// A structural deletion of 25000 to 26000 lies 5000 bases upstream of the
    /// reverse transcript ending at 20000: structural rows are annotated, the
    /// batch dispatch handing the plugin the structural variant because
    /// TSSDistance.pm:73-75 declares `BaseVariationFeature`, where a plugin
    /// with Ensembl VEP's default `variant_feature_types` would be skipped
    /// (`OutputFactory.pm:2229`).
    #[test]
    fn a_structural_variant_row_is_annotated() {
        let plugin = plugin(&[]);
        let mut deletion = span(25000, 26000);
        deletion.ref_allele = b"N".to_vec();
        deletion.alt_alleles = vec![b"<DEL>".to_vec()];
        deletion.is_structural = true;
        deletion.sv_end = Some(26000);
        assert_eq!(
            distance(&plugin, &row(-1, 10000, 20000), &deletion).as_deref(),
            Some("5000")
        );

        assert!(plugin.annotates_structural_variants());
        deletion.transcript_consequences.push(row(-1, 10000, 20000));
        plugin
            .run_batch(std::slice::from_mut(&mut deletion))
            .unwrap();
        assert_eq!(
            deletion.transcript_consequences[0]
                .plugin_data
                .get(OUTPUT_FIELD)
                .map(String::as_str),
            Some("5000")
        );
    }

    /// `params_to_hash` reads the value by its truth: `0` and the empty value are
    /// off, any other value on.
    #[test]
    fn both_direction_follows_perl_truth() {
        assert!(!plugin(&[]).both_direction);
        assert!(!plugin(&["both_direction=0"]).both_direction);
        assert!(!plugin(&["both_direction="]).both_direction);
        assert!(plugin(&["both_direction=1"]).both_direction);
        assert!(plugin(&["both_direction=yes"]).both_direction);
        assert!(plugin(&["both_direction=00"]).both_direction);
    }

    /// A token without `=` leaves `params_to_hash` empty (BaseVepPlugin.pm:123),
    /// so Ensembl VEP's plugin runs with `both_direction` off; the token is not
    /// read here either. A key the plugin does not declare is refused.
    #[test]
    fn a_bare_token_is_not_read_and_an_unknown_key_is_refused() {
        assert!(!plugin(&["1"]).both_direction);
        assert!(!plugin(&["both_direction"]).both_direction);
        let mut plugin = TssDistancePlugin::new();
        let error = plugin
            .init(&["direction=1".to_string()])
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            "plugin initialization failed: TSSDistance: unknown parameter 'direction'; \
             TSSDistance accepts both_direction"
        );
    }
}
