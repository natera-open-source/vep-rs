// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! NMD plugin: the nonsense-mediated mRNA decay escape prediction for a
//! truncating variant, as Ensembl VEP's `NMD.pm` makes it.
//!
//! **Usage:** `--plugin NMD` (the plugin takes no parameter; a positional one is
//! not read, and is logged; a `key=value` one is refused)
//!
//! Output field: `NMD`, with the value `NMD_escaping_variant` on a transcript
//! row whose terms include `stop_gained`, `frameshift_variant`,
//! `splice_donor_variant` or `splice_acceptor_variant` (`NMD.pm:65, :83`),
//! whose CDS position has both ends defined (`:103-105`) and which meets one of
//! four rules: the variant ends inside the transcript's last exon
//! (`:142-150`); it ends inside the 51 bases before the 3' end of the
//! penultimate exon (`:158-169`); its CDS end is at or before coding base 101
//! (`:98`); or the transcript has one exon (`:109`). Every other row carries
//! nothing, a structural variant's rows included: Ensembl VEP hands a plugin
//! only the variant classes its `variant_feature_types` names
//! (`OutputFactory.pm:2220`), and `NMD.pm` keeps the default, sequence
//! variants alone (`BaseVepPlugin.pm:100`).
//!
//! Exon coordinates come from [`TranscriptConsequence::loftee_ctx`], the
//! structural context the pipeline attaches to every row when this plugin is
//! loaded; the CDS end is read back from the row's `CDS_position` column, or
//! for an allele of a `--minimal` split from the allele's own span
//! ([`TranscriptConsequence::own_span_positions`]), the span its genomic end
//! and terms describe.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`]
//! via the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use indexmap::IndexMap;

use crate::params::parse_none;
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::{Consequence, LofteeContext, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "NMD";

/// The one value the field takes (`NMD.pm:67`).
const ESCAPING: &str = "NMD_escaping_variant";

/// The terms a row needs one of to be evaluated (`NMD.pm:65`).
const GATE_TERMS: [Consequence; 4] = [
    Consequence::StopGained,
    Consequence::FrameshiftVariant,
    Consequence::SpliceDonorVariant,
    Consequence::SpliceAcceptorVariant,
];

/// A variant whose CDS end is at or before this coding base escapes
/// (`NMD.pm:98`).
const FIRST_CODING_BASES: u64 = 101;

/// Distance from the penultimate exon's 3' end, in bases, inside which the
/// variant's end escapes (`NMD.pm:160, :166`).
const PENULTIMATE_EXON_WINDOW: u64 = 51;

/// Nonsense-mediated decay escape prediction plugin.
#[derive(Debug, Default)]
pub struct NmdPlugin;

impl NmdPlugin {
    /// Create the plugin; `--plugin NMD` carries no parameter to initialize.
    pub fn new() -> Self {
        Self
    }
}

impl BuiltinPlugin for NmdPlugin {
    fn name(&self) -> &str {
        "NMD"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "Nonsense-mediated mRNA decay escaping variants prediction".to_string(),
        )]
    }

    /// A positional parameter is not read, as the Ensembl VEP module reads
    /// none (NMD.pm defines no `new`); it is logged. A `key=value` list is
    /// refused as every built-in refuses an undeclared key.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        parse_none(self.name(), params)
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if escapes_nmd(consequence, variant) {
            result.insert(OUTPUT_FIELD.to_string(), ESCAPING.to_string());
        }
        Ok(result)
    }
}

/// `NMD.pm` `run` (`:79-117`): a gate term, both CDS ends defined, and one of
/// the four rules, on a sequence variant's row alone.
fn escapes_nmd(consequence: &TranscriptConsequence, variant: &InputVariant) -> bool {
    // A structural deletion inside one coding exon carries `frameshift_variant`
    // and a CDS position too, but Ensembl VEP never runs this plugin on a
    // `StructuralVariationFeature` (OutputFactory.pm:2220, BaseVepPlugin.pm:100).
    if variant.is_structural {
        return false;
    }
    if !consequence
        .consequences
        .iter()
        .any(|term| GATE_TERMS.contains(term))
    {
        return false;
    }
    let Some(cds_end) = cds_end(consequence, variant) else {
        return false;
    };
    let Some(ctx) = consequence.loftee_ctx.as_deref() else {
        return false;
    };
    cds_end <= FIRST_CODING_BASES
        || ctx.exons.len() == 1
        || ends_in_last_or_penultimate_exon(ctx, variant.end)
}

/// The row's CDS end as `TranscriptVariation::cds_end` reads it
/// (`BaseTranscriptVariation.pm:252-294`): `None` unless both ends of the
/// `CDS_position` column are coordinates (`?` marks a gap). The column prints
/// an insertion's two coordinates swapped (`Utils::format_coords`), so its CDS
/// end is the first number. For an allele of a `--minimal` split the column
/// prints the record's span and the allele's own span is `own_span_positions`,
/// the one judged, as `variant.end` is the allele's own end: an own span with
/// no CDS coordinate is `None`, whatever the record's column prints.
fn cds_end(consequence: &TranscriptConsequence, variant: &InputVariant) -> Option<u64> {
    let column = match consequence.own_span_positions.as_deref() {
        Some(own) => own.cds.as_deref()?,
        None => consequence.cds_position.as_deref()?,
    };
    let (first, last) = column.split_once('-').unwrap_or((column, column));
    let first: u64 = first.parse().ok()?;
    let last: u64 = last.parse().ok()?;
    let is_insertion = variant.end < variant.start;
    Some(if is_insertion { first } else { last })
}

/// `variant_exon_check` (`NMD.pm:132-176`): the variant's end lies in the
/// transcript's last exon, or in its penultimate exon within
/// [`PENULTIMATE_EXON_WINDOW`] bases of that exon's 3' end, the exon's end on
/// the forward strand and its start on the reverse. Exons are in transcript
/// order.
fn ends_in_last_or_penultimate_exon(ctx: &LofteeContext, variant_end: u64) -> bool {
    let Some(last) = ctx.exons.last() else {
        return false;
    };
    if (last.start..=last.end).contains(&variant_end) {
        return true;
    }
    let Some(penultimate) = ctx.exons.iter().rev().nth(1) else {
        return false;
    };
    let window = if ctx.strand == -1 {
        penultimate.start..=penultimate.start + PENULTIMATE_EXON_WINDOW
    } else {
        penultimate.end.saturating_sub(PENULTIMATE_EXON_WINDOW)..=penultimate.end
    };
    window.contains(&variant_end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::{FeatureType, OwnSpanPositions};
    use vep_core::transcript::Exon;
    use vep_core::variant::VariantClass;

    /// A transcript of `exons`, given as genomic `(start, end)` pairs in
    /// transcript order, on `strand`.
    fn context(strand: i8, exons: &[(u64, u64)]) -> LofteeContext {
        LofteeContext {
            exons: exons
                .iter()
                .enumerate()
                .map(|(i, &(start, end))| Exon {
                    start,
                    end,
                    rank: i as u32 + 1,
                    phase: 0,
                    end_phase: 0,
                })
                .collect(),
            introns: Vec::new(),
            strand,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            cds_start_nf: false,
            cds_end_nf: false,
        }
    }

    /// A transcript row carrying `terms` and the `CDS_position` column `cds`.
    fn row(terms: &[Consequence], cds: Option<&str>, ctx: LofteeContext) -> TranscriptConsequence {
        TranscriptConsequence {
            consequences: terms.iter().copied().collect(),
            cds_position: cds.map(str::to_string),
            strand: ctx.strand,
            feature_type: FeatureType::Transcript,
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        }
    }

    /// A single-base substitution at `pos`.
    fn snv(pos: u64) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, b"C".to_vec(), b"T".to_vec())
    }

    fn nmd(consequence: &TranscriptConsequence, variant: &InputVariant) -> Option<String> {
        NmdPlugin::new()
            .run(consequence, variant)
            .unwrap()
            .get(OUTPUT_FIELD)
            .cloned()
    }

    /// Three exons on the forward strand, in transcript order.
    const FORWARD_THREE: [(u64, u64); 3] = [(1000, 1100), (2000, 2100), (3000, 3100)];

    /// The same three exons on the reverse strand, where transcript order runs
    /// from the rightmost exon to the leftmost.
    const REVERSE_THREE: [(u64, u64); 3] = [(3000, 3100), (2000, 2100), (1000, 1100)];

    #[test]
    fn only_the_four_gate_terms_are_evaluated() {
        // A stop in the last exon with a CDS end past 101: the last-exon rule
        // alone would fire, so every outcome below is the gate's.
        let last_exon = snv(3050);
        for term in GATE_TERMS {
            let tc = row(&[term], Some("500"), context(1, &FORWARD_THREE));
            assert_eq!(nmd(&tc, &last_exon).as_deref(), Some(ESCAPING), "{term:?}");
        }
        let tc = row(
            &[Consequence::MissenseVariant, Consequence::StopGained],
            Some("500"),
            context(1, &FORWARD_THREE),
        );
        assert_eq!(nmd(&tc, &last_exon).as_deref(), Some(ESCAPING));
        for term in [
            Consequence::MissenseVariant,
            Consequence::StopLost,
            Consequence::SpliceRegionVariant,
            Consequence::InframeDeletion,
            Consequence::StartLost,
        ] {
            let tc = row(&[term], Some("500"), context(1, &FORWARD_THREE));
            assert_eq!(nmd(&tc, &last_exon), None, "{term:?}");
        }
    }

    #[test]
    fn no_cds_coordinates_means_no_annotation() {
        // The same last-exon stop with the CDS column absent, open at either
        // end (a gap on one side of the mapping), and present.
        let last_exon = snv(3050);
        for cds in [None, Some("?-500"), Some("500-?")] {
            let tc = row(&[Consequence::StopGained], cds, context(1, &FORWARD_THREE));
            assert_eq!(nmd(&tc, &last_exon), None, "{cds:?}");
        }
        let tc = row(
            &[Consequence::StopGained],
            Some("500-502"),
            context(1, &FORWARD_THREE),
        );
        assert_eq!(nmd(&tc, &last_exon).as_deref(), Some(ESCAPING));
    }

    #[test]
    fn cds_end_at_101_escapes_and_102_does_not() {
        // A stop in the first of three exons, 1,050 bases before the
        // penultimate exon's end: only the first-coding-bases rule can fire.
        let first_exon = snv(1050);
        let ctx = || context(1, &FORWARD_THREE);
        assert_eq!(
            nmd(
                &row(&[Consequence::StopGained], Some("101"), ctx()),
                &first_exon
            )
            .as_deref(),
            Some(ESCAPING)
        );
        assert_eq!(
            nmd(
                &row(&[Consequence::StopGained], Some("102"), ctx()),
                &first_exon
            ),
            None
        );
        // A range's end is its second number (format_coords prints start-end).
        let deletion = InputVariant::new("21".into(), 1049, 1050, b"AT".to_vec(), b"-".to_vec());
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], Some("99-101"), ctx()),
                &deletion
            )
            .as_deref(),
            Some(ESCAPING)
        );
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], Some("100-102"), ctx()),
                &deletion
            ),
            None
        );
        // An insertion's column is printed end-start, so its CDS end is the
        // first number: 101-102 is an insertion after coding base 101.
        let insertion = InputVariant::new("21".into(), 1051, 1050, b"-".to_vec(), b"T".to_vec());
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], Some("101-102"), ctx()),
                &insertion
            )
            .as_deref(),
            Some(ESCAPING)
        );
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], Some("102-103"), ctx()),
                &insertion
            ),
            None
        );
    }

    #[test]
    fn last_exon_hit_escapes_on_both_strands() {
        // CDS end 500 keeps the first-coding-bases rule out; the variant sits
        // more than 51 bases from every penultimate-exon junction.
        let cds = Some("500");
        // Forward: the last exon in transcript order is the rightmost.
        let forward = || context(1, &FORWARD_THREE);
        for pos in [3000, 3050, 3100] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, forward()), &snv(pos)).as_deref(),
                Some(ESCAPING),
                "forward {pos}"
            );
        }
        for pos in [2999, 3101, 1050] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, forward()), &snv(pos)),
                None,
                "forward {pos}"
            );
        }
        // Reverse: the last exon in transcript order is the leftmost.
        let reverse = || context(-1, &REVERSE_THREE);
        for pos in [1000, 1050, 1100] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, reverse()), &snv(pos)).as_deref(),
                Some(ESCAPING),
                "reverse {pos}"
            );
        }
        for pos in [999, 1101, 3050] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, reverse()), &snv(pos)),
                None,
                "reverse {pos}"
            );
        }
        // A deletion is placed by its end: one ending in the last exon escapes,
        // one ending just before it does not.
        let into_last = InputVariant::new("21".into(), 2990, 3002, b"A".repeat(13), b"-".to_vec());
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], cds, forward()),
                &into_last
            )
            .as_deref(),
            Some(ESCAPING)
        );
        let before_last =
            InputVariant::new("21".into(), 2990, 2999, b"A".repeat(10), b"-".to_vec());
        assert_eq!(
            nmd(
                &row(&[Consequence::FrameshiftVariant], cds, forward()),
                &before_last
            ),
            None
        );
    }

    #[test]
    fn penultimate_exon_window_is_51_bases_from_the_junction_on_both_strands() {
        let cds = Some("500");
        // Forward: the window is the penultimate exon's last 52 bases,
        // 2049-2100.
        let forward = || context(1, &FORWARD_THREE);
        for pos in [2049, 2075, 2100] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, forward()), &snv(pos)).as_deref(),
                Some(ESCAPING),
                "forward {pos}"
            );
        }
        for pos in [2048, 2000, 1100] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, forward()), &snv(pos)),
                None,
                "forward {pos}"
            );
        }
        // Reverse: the penultimate exon is 2000-2100 and its 3' end is its
        // start, so the window is 2000-2051.
        let reverse = || context(-1, &REVERSE_THREE);
        for pos in [2000, 2025, 2051] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, reverse()), &snv(pos)).as_deref(),
                Some(ESCAPING),
                "reverse {pos}"
            );
        }
        for pos in [2052, 2100, 3050] {
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, reverse()), &snv(pos)),
                None,
                "reverse {pos}"
            );
        }
        // The window measures from the junction of the penultimate exon alone:
        // the same offset in the first of four exons does not escape.
        let four = context(1, &[(1000, 1100), (2000, 2100), (3000, 3100), (4000, 4100)]);
        assert_eq!(
            nmd(&row(&[Consequence::StopGained], cds, four), &snv(2075)),
            None
        );
        // The penultimate window is still read when the window reaches past
        // the exon's start: exon 1 is 50 bases, its whole span escapes.
        let short = context(1, &[(1000, 1049), (2000, 2100)]);
        assert_eq!(
            nmd(&row(&[Consequence::StopGained], cds, short), &snv(1000)).as_deref(),
            Some(ESCAPING)
        );
    }

    #[test]
    fn intronless_transcript_escapes() {
        let cds = Some("500");
        for strand in [1, -1] {
            let one = || context(strand, &[(1000, 2000)]);
            assert_eq!(
                nmd(&row(&[Consequence::StopGained], cds, one()), &snv(1500)).as_deref(),
                Some(ESCAPING),
                "strand {strand}"
            );
            assert_eq!(
                nmd(
                    &row(&[Consequence::FrameshiftVariant], cds, one()),
                    &snv(1999)
                )
                .as_deref(),
                Some(ESCAPING),
                "strand {strand}"
            );
        }
        // Two exons with the stop in the first, outside the penultimate
        // window: no rule fires.
        let two = context(1, &[(1000, 2000), (3000, 3100)]);
        assert_eq!(
            nmd(&row(&[Consequence::StopGained], cds, two), &snv(1500)),
            None
        );
    }

    /// The `NMD` key is plugin data on the row: it joins the `Extra` column and
    /// leaves the row's consequence terms as they were (`SpliceRegion.pm`'s
    /// shape, where a plugin adds a term, is not this plugin's).
    #[test]
    fn nmd_is_an_extra_key_not_a_consequence_term() {
        let mut variant = snv(3050);
        variant.transcript_consequences.push(row(
            &[Consequence::StopGained],
            Some("500"),
            context(1, &FORWARD_THREE),
        ));
        variant.transcript_consequences.push(row(
            &[Consequence::MissenseVariant],
            Some("500"),
            context(1, &FORWARD_THREE),
        ));
        NmdPlugin::new()
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();

        let [escaping, missense] = &variant.transcript_consequences[..] else {
            panic!("two rows");
        };
        assert_eq!(
            escaping.plugin_data.get(OUTPUT_FIELD).map(String::as_str),
            Some(ESCAPING)
        );
        assert_eq!(
            escaping.consequences.as_slice(),
            [Consequence::StopGained].as_slice()
        );
        assert!(missense.plugin_data.is_empty());
        assert_eq!(
            missense.consequences.as_slice(),
            [Consequence::MissenseVariant].as_slice()
        );
        assert!(variant.plugin_data.is_empty());
        assert_eq!(
            NmdPlugin::new().header_info(),
            [(
                OUTPUT_FIELD.to_string(),
                "Nonsense-mediated mRNA decay escaping variants prediction".to_string()
            )]
        );
    }

    /// The Ensembl VEP module reads no parameter (NMD.pm defines no `new`, so
    /// `BaseVepPlugin::new` stores the list unread), and a positional one is
    /// not read here either; a `key=value` list is refused as the undeclared
    /// key every built-in refuses.
    #[test]
    fn a_parameter_is_not_read() {
        let mut plugin = NmdPlugin::new();
        assert!(plugin.init(&[]).is_ok());
        assert!(plugin.init(&["1".to_string()]).is_ok());
        assert!(plugin.init(&["a".to_string(), "b".to_string()]).is_ok());
        assert_eq!(
            plugin
                .init(&["file=x".to_string()])
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: NMD: unknown parameter 'file'; NMD takes no parameter"
        );
    }

    /// A structural deletion inside the last exon, `frameshift_variant` and
    /// CDS position included, carries nothing; the same span as a sequence
    /// deletion escapes by the last-exon rule.
    #[test]
    fn a_structural_variant_row_carries_nothing() {
        let last_exon_frameshift = || {
            row(
                &[Consequence::FrameshiftVariant],
                Some("413-416"),
                context(1, &FORWARD_THREE),
            )
        };
        let mut structural =
            InputVariant::new("21".into(), 3001, 3004, b"N".to_vec(), b"<DEL>".to_vec());
        structural.variant_class = VariantClass::StructuralDeletion;
        structural.is_structural = true;
        structural.sv_end = Some(3004);
        assert_eq!(nmd(&last_exon_frameshift(), &structural), None);

        let sequence = InputVariant::new("21".into(), 3001, 3004, b"ACGT".to_vec(), b"-".to_vec());
        assert_eq!(
            nmd(&last_exon_frameshift(), &sequence).as_deref(),
            Some(ESCAPING)
        );
    }

    /// A row without the structural context carries nothing: the plugin never
    /// guesses from the summary fields.
    #[test]
    fn a_row_without_exon_structure_carries_nothing() {
        let mut tc = row(
            &[Consequence::StopGained],
            Some("50"),
            context(1, &FORWARD_THREE),
        );
        tc.loftee_ctx = None;
        assert_eq!(nmd(&tc, &snv(1050)), None);
    }

    /// An allele of a `--minimal` split is judged over its own trimmed span:
    /// the `CDS_position` column prints the record's span and the allele's own
    /// is `own_span_positions`. The record `GTG GG,G` splits into a deletion of
    /// coding base 101 (the record's span is 101-102) and a deletion of 101-102.
    #[test]
    fn a_minimal_split_allele_is_judged_over_its_own_cds_span() {
        let record_row = || {
            row(
                &[Consequence::FrameshiftVariant],
                Some("101-102"),
                context(1, &FORWARD_THREE),
            )
        };
        let with_own_span = |mut tc: TranscriptConsequence, cds: Option<&str>| {
            tc.own_span_positions = Some(Box::new(OwnSpanPositions {
                cdna: None,
                cds: cds.map(str::to_string),
                protein: None,
            }));
            tc
        };
        let mut one_base = InputVariant::new("21".into(), 1050, 1050, b"T".to_vec(), b"-".to_vec());
        one_base.original_start = Some(1050);
        one_base.original_end = Some(1051);
        assert_eq!(
            nmd(&with_own_span(record_row(), Some("101")), &one_base).as_deref(),
            Some(ESCAPING)
        );
        let mut two_bases =
            InputVariant::new("21".into(), 1050, 1051, b"TG".to_vec(), b"-".to_vec());
        two_bases.original_start = Some(1050);
        two_bases.original_end = Some(1051);
        assert_eq!(
            nmd(&with_own_span(record_row(), Some("101-102")), &two_bases),
            None
        );
        // The column stands where the allele's own span is the record's.
        assert_eq!(nmd(&record_row(), &two_bases), None);
        // An own span with no CDS coordinate at either end (a whole exon deleted
        // from intron to intron) carries nothing, however the record's column,
        // here a span the first rule would mark, reads.
        let whole_exon = with_own_span(
            row(
                &[
                    Consequence::SpliceAcceptorVariant,
                    Consequence::SpliceDonorVariant,
                ],
                Some("101"),
                context(1, &FORWARD_THREE),
            ),
            None,
        );
        let mut intron_to_intron =
            InputVariant::new("21".into(), 1950, 2150, b"N".to_vec(), b"-".to_vec());
        intron_to_intron.original_start = Some(1050);
        intron_to_intron.original_end = Some(3050);
        assert_eq!(nmd(&whole_exon, &intron_to_intron), None);
    }
}
