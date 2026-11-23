// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! NearestExonJB plugin: the nearest exon junction boundary of an exonic
//! variant, as Ensembl VEP's `NearestExonJB.pm` writes it.
//!
//! **Usage:** `--plugin NearestExonJB[,max_range=<bases>][,intronic=1]`
//!
//! Output field: `NearestExonJB`, one `ExonID|distance|start/end|length` item
//! per reported exon, items joined by `,` and their parts by `|`, or by `+`
//! under the VCF output format (`NearestExonJB.pm:80, :98`). For each exon the
//! variant overlaps (`:184`; the transcript's exons stretched by 12 bases when
//! it has an intron of at most 13 bases, as `_overlapped_exons` stretches them)
//! the plugin measures the variant's start to the exon's start and its end to
//! the exon's end, keeps the shorter as the exon's distance with the side it
//! lies on, `start_end` on a tie (`:225-239`), and reports the exons whose
//! distance is the smallest, when that is within `max_range` (default 10000,
//! `:83, :186, :267-271`). Under `intronic=1` a row with no overlapped exon
//! whose terms include `intron_variant` reports the two exons flanking the
//! intron, each with its own distance, `max_range` not applied (`:190-217,
//! :244-264`); a run of two overlapped exons reports both the same way (`:244`).
//! A row overlapping no exon carries nothing, as does a structural variant's
//! row (`BaseVepPlugin.pm:100`). Ensembl VEP lists several items in hash order;
//! vep-rs lists them in transcript order.
//!
//! Both values are read as Perl's numeric comparisons read their text:
//! `$CONFIG{intronic} == 1` (`:135, :190, :244`) and `$min > $startD`
//! (`:230-238`) take the leading number, 0 when none leads it, so `intronic=01`
//! is on, `intronic=2` and `intronic=` are off, `max_range=50.0` admits a
//! distance of 49, and `max_range=abc` admits none. When no exon lies nearer
//! than `max_range`, `$min` still holds the value's text and
//! `exists $dists{$exon}{$min}` (`:268`) matches it against the distances'
//! decimal text, so an exon exactly at the range is reported under
//! `max_range=50` and the default, and not under `50.0`, `050` or `+50`. The
//! same hash keys each exon's length under `len` (`:229-237`), so `max_range=len`
//! makes Ensembl VEP report every overlapped exon with `len` as its distance;
//! vep-rs reads `len` as 0, like any other non-numeric value, and reports
//! nothing, the documented class `nearest_exon_jb_len_key_collision`.
//!
//! The exon model comes from [`TranscriptConsequence::loftee_ctx`], which the
//! pipeline attaches to every row when this plugin is loaded, the exon ids from
//! [`TranscriptConsequence::exon_ids`], which the cache loads for this plugin
//! alone, and the intron number from the row's `INTRON` column. An exon without
//! a stable id (a RefSeq aligned transcript's, written under `--all_refseq`) is
//! reported with the empty id, `|<distance>|<side>|<length>`, as Ensembl VEP
//! writes it (`:220, :269`); where Ensembl VEP keys its exons by id and so folds
//! every id-less exon of a row onto one item (`:228-237`), vep-rs reports each
//! exon's own item, the documented class `nearest_exon_jb_id_less_exon_fold`.
//!
//! `gff3=<file>`, which reads exon locations from a GFF3 file (`:52-58`), is
//! refused: vep-rs reads the cache's transcript model. A `gff3` value Perl
//! reads as false, `0` or the empty string, is no file (`:115`), and
//! `parent_only`, `0` or `1` (`:106-108`), is read by a GFF3 file alone, so
//! both leave the cache's model in use.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`]
//! via the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use std::cmp::Ordering;
use std::fmt::Write;

use indexmap::IndexMap;

use crate::params::{perl_numeric, PluginParams};
use crate::traits::{BuiltinPlugin, OutputFormat, PluginError};
use vep_core::consequence::{Consequence, LofteeContext, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::transcript::Exon;
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "NearestExonJB";

/// The parameters the plugin reads (`NearestExonJB.pm:46-58`).
const PARAMETERS: [&str; 4] = ["max_range", "intronic", "gff3", "parent_only"];

/// The default search range in bases (`NearestExonJB.pm:83`).
const DEFAULT_MAX_RANGE: u64 = 10000;

/// Bases each exon stretches by on a transcript with a frameshift intron
/// (`BaseTranscriptVariation.pm:861`).
const FRAMESHIFT_INTRON_STRETCH: u64 = 12;

/// Largest `end - start` of an intron that counts as a frameshift intron, one
/// of at most 13 bases (`BaseTranscriptVariation.pm:1032`); an intron whose
/// start is past its end, as abutting or overlapping exons leave it, is skipped
/// before the test (`:1014`).
const FRAMESHIFT_INTRON_SPAN: u64 = 12;

/// Nearest exon junction boundary plugin.
#[derive(Debug)]
pub struct NearestExonJbPlugin {
    /// `max_range` as `$min > $startD` numifies its text (`NearestExonJB.pm:186,
    /// :230-238`).
    max_range: f64,
    /// The distance reported when no exon lies nearer than `max_range`: `$min`
    /// then still holds the value's text, which `exists $dists{$exon}{$min}`
    /// (`:268`) matches against the distances' decimal text, so only a value
    /// written as a plain decimal integer (`50`, not `50.0`, `050` or `+50`)
    /// reports an exon exactly at the range. The default 10000 does.
    at_range: Option<u64>,
    /// `intronic` as `$CONFIG{intronic} == 1` numifies its text
    /// (`NearestExonJB.pm:135, :190, :244`): on when that number is 1.
    intronic: bool,
    separator: char,
}

impl Default for NearestExonJbPlugin {
    fn default() -> Self {
        Self {
            max_range: DEFAULT_MAX_RANGE as f64,
            at_range: Some(DEFAULT_MAX_RANGE),
            intronic: false,
            separator: '|',
        }
    }
}

impl NearestExonJbPlugin {
    /// Create the plugin; `init` reads `max_range` and `intronic`.
    pub fn new() -> Self {
        Self::default()
    }
}

impl crate::traits::sealed::Sealed for NearestExonJbPlugin {}

impl BuiltinPlugin for NearestExonJbPlugin {
    fn name(&self) -> &str {
        "NearestExonJB"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        let scope = if self.intronic {
            ""
        } else {
            " (coding sequence variants only)"
        };
        let parts = ["ExonID", "distance", "start/end", "length"].join(&self.separator.to_string());
        vec![(
            OUTPUT_FIELD.to_string(),
            format!("Nearest Exon Junction Boundary{scope}. Format:{parts}"),
        )]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let parsed = PluginParams::parse(self.name(), &PARAMETERS, params)?;
        if let Some(token) = parsed.positional().first() {
            return Err(PluginError::Init(format!(
                "NearestExonJB: parameter '{token}' is not key=value; \
                 use max_range=<bases> or intronic=1"
            )));
        }
        // `if($self->{gff3_file})` (`NearestExonJB.pm:115`): a value Perl reads
        // as false, `0` or the empty string, is no file.
        if parsed.get("gff3").is_some_and(perl_true) {
            return Err(PluginError::Init(
                "NearestExonJB: gff3= reads exon locations from a GFF3 file; vep-rs reads \
                 them from the cache's transcript model and takes no file"
                    .to_string(),
            ));
        }
        // `parent_only` is read and checked whether or not a GFF3 file is given
        // (`:106-108`); without one nothing reads it.
        if let Some(value) = parsed.get("parent_only") {
            if !matches!(value, "0" | "1") {
                return Err(PluginError::Init(format!(
                    "NearestExonJB: parent_only must be 0 or 1, got '{value}'"
                )));
            }
        }
        if let Some(value) = parsed.get("max_range") {
            self.max_range = perl_numeric(self.name(), "max_range", value);
            self.at_range = value
                .parse()
                .ok()
                .filter(|range: &u64| range.to_string() == value);
        }
        if let Some(value) = parsed.get("intronic") {
            self.intronic = perl_numeric(self.name(), "intronic", value) == 1.0;
        }
        Ok(())
    }

    fn set_output_format(&mut self, format: OutputFormat) {
        self.separator = if format == OutputFormat::Vcf {
            '+'
        } else {
            '|'
        };
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if variant.is_structural {
            return Ok(result);
        }
        let Some(ctx) = consequence.loftee_ctx.as_deref() else {
            return Ok(result);
        };
        let items = self.items(ctx, consequence, variant);
        if !items.is_empty() {
            result.insert(OUTPUT_FIELD.to_string(), items);
        }
        Ok(result)
    }
}

/// An exon's junction nearest the variant (`NearestExonJB.pm:225-239`); `index`
/// is the exon's position in the transcript model, the index of its id in
/// [`TranscriptConsequence::exon_ids`].
struct NearestJunction<'a> {
    exon: &'a Exon,
    index: usize,
    distance: u64,
    side: &'static str,
}

impl NearestExonJbPlugin {
    /// `NearestExonJB.pm` `run` (`:173-277`): the reported items in transcript
    /// order, joined by `,` (`:274`); empty when none is reported.
    fn items(
        &self,
        ctx: &LofteeContext,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> String {
        let mut exons = overlapped_exons(ctx, variant);
        if exons.is_empty() && self.intronic {
            exons = flanking_exons(ctx, consequence);
        }
        let junctions: Vec<NearestJunction<'_>> = exons
            .into_iter()
            .map(|(index, exon)| nearest_junction(exon, index, variant))
            .collect();
        // Two exons under `intronic` are both reported; otherwise `$min` starts
        // at `max_range`, each exon's distance below it lowers it, and the exons
        // at it are reported (`NearestExonJB.pm:186, :230-238, :268`): the
        // nearest when one lies within the range, else the exons exactly at a
        // range whose text is a plain decimal integer.
        let report_both = junctions.len() == 2 && self.intronic;
        let nearest = junctions.iter().map(|junction| junction.distance).min();
        let reported = match nearest {
            Some(distance) if (distance as f64) < self.max_range => Some(distance),
            _ => self.at_range,
        };
        let ids = consequence.exon_ids.as_deref().unwrap_or(&[]);
        let mut items = String::new();
        for junction in junctions
            .iter()
            .filter(|junction| report_both || Some(junction.distance) == reported)
        {
            if !items.is_empty() {
                items.push(',');
            }
            // An exon without an id keeps its item, written with the empty id
            // (`NearestExonJB.pm:220, :269`).
            let id = ids.get(junction.index).map_or("", |id| &**id);
            let length = junction.exon.end - junction.exon.start + 1;
            let (sep, distance, side) = (self.separator, junction.distance, junction.side);
            write!(items, "{id}{sep}{distance}{sep}{side}{sep}{length}")
                .expect("writing to a String");
        }
        items
    }
}

/// The span the exon set is read over and the distances measure from
/// (`_overlapped_exons`, `NearestExonJB.pm:225-226`): the variant's own. An
/// allele of a `--minimal` split record reads its own trimmed span (the
/// documented classes `minimal_plugin_sees_the_rejoined_record` and, for the
/// item string Ensembl VEP's cache keyed by transcript and span hands the
/// record's other rows, `NearestExonJB.pm:179-181`,
/// `nearest_exon_jb_shared_span_cache`). The row's `shift_length` moves its
/// terms, never this span, where Ensembl VEP's plugin reads the memo its
/// consequence tests filled at the shifted span, the documented class
/// `nearest_exon_jb_shifted_span`.
fn plugin_span(variant: &InputVariant) -> (u64, u64) {
    (variant.start, variant.end)
}

/// `_overlapped_exons` (`BaseTranscriptVariation.pm:848-870`): the exons
/// overlapping the span, each stretched by 12 bases when the transcript has an
/// intron of at most 13 bases.
fn overlapped_exons<'a>(ctx: &'a LofteeContext, variant: &InputVariant) -> Vec<(usize, &'a Exon)> {
    let stretch = if ctx.introns.iter().any(|intron| {
        intron.start <= intron.end && intron.end - intron.start <= FRAMESHIFT_INTRON_SPAN
    }) {
        FRAMESHIFT_INTRON_STRETCH
    } else {
        0
    };
    let (start, end) = plugin_span(variant);
    let (low, high) = (start.min(end), start.max(end));
    ctx.exons
        .iter()
        .enumerate()
        .filter(|(_, exon)| exon.start <= high + stretch && low <= exon.end + stretch)
        .collect()
}

/// The two exons flanking the intron an `intron_variant` row lies in
/// (`NearestExonJB.pm:190-217`): the `INTRON` column's first number names the
/// intron in transcript order, so exons `n` and `n + 1` of the model flank it;
/// empty for any other row.
fn flanking_exons<'a>(
    ctx: &'a LofteeContext,
    consequence: &TranscriptConsequence,
) -> Vec<(usize, &'a Exon)> {
    let Some(intron) = consequence.intron.as_deref() else {
        return Vec::new();
    };
    if !consequence
        .consequences
        .contains(&Consequence::IntronVariant)
    {
        return Vec::new();
    }
    // Perl numifies the column's first field, so `3-4/10` reads as intron 3.
    let digits: String = intron.chars().take_while(char::is_ascii_digit).collect();
    let Ok(number) = digits.parse::<usize>() else {
        return Vec::new();
    };
    let Some(before) = number.checked_sub(1) else {
        return Vec::new();
    };
    match (ctx.exons.get(before), ctx.exons.get(number)) {
        (Some(exon_before), Some(exon_after)) => {
            vec![(before, exon_before), (number, exon_after)]
        }
        _ => Vec::new(),
    }
}

/// `NearestExonJB.pm:225-239`: the shorter of the span's start to the exon's
/// start and its end to the exon's end, with its side.
fn nearest_junction<'a>(
    exon: &'a Exon,
    index: usize,
    variant: &InputVariant,
) -> NearestJunction<'a> {
    let (start, end) = plugin_span(variant);
    let start_distance = start.abs_diff(exon.start);
    let end_distance = end.abs_diff(exon.end);
    let (distance, side) = match start_distance.cmp(&end_distance) {
        Ordering::Less => (start_distance, "start"),
        Ordering::Greater => (end_distance, "end"),
        Ordering::Equal => (start_distance, "start_end"),
    };
    NearestJunction {
        exon,
        index,
        distance,
        side,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Intron;

    /// A transcript of `exons`, given as genomic `(start, end)` pairs in
    /// transcript order, on `strand`, with the introns between them.
    fn context(strand: i8, exons: &[(u64, u64)]) -> LofteeContext {
        let exons: Vec<Exon> = exons
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| Exon {
                start,
                end,
                rank: i as u32 + 1,
                phase: 0,
                end_phase: 0,
            })
            .collect();
        let introns = exons
            .windows(2)
            .enumerate()
            .map(|(i, pair)| {
                let (a, b) = (&pair[0], &pair[1]);
                let (left, right) = if a.start < b.start { (a, b) } else { (b, a) };
                Intron {
                    start: left.end + 1,
                    end: right.start - 1,
                    rank: i as u32 + 1,
                }
            })
            .collect();
        let mut context =
            LofteeContext::from_transcript(&vep_core::transcript::Transcript::default());
        context.exons = exons;
        context.introns = introns;
        context.strand = strand;
        context.coding_region_start = None;
        context.coding_region_end = None;
        context.translation_start = None;
        context.translation_end = None;
        context.cdna_coding_start = None;
        context.cdna_coding_end = None;
        context.cds_start_nf = false;
        context.cds_end_nf = false;
        context
    }

    fn row(
        terms: &[Consequence],
        intron: Option<&str>,
        ctx: LofteeContext,
    ) -> TranscriptConsequence {
        // The ids the cache would carry: `ENSE` plus the exon's rank, padded.
        let exon_ids: Arc<[Arc<str>]> = ctx
            .exons
            .iter()
            .map(|exon| Arc::from(format!("ENSE{:011}", exon.rank)))
            .collect();
        let mut consequence = TranscriptConsequence::default();
        consequence.consequences = terms.iter().copied().collect();
        consequence.intron = intron.map(str::to_string);
        consequence.strand = ctx.strand;
        consequence.feature_type = FeatureType::Transcript;
        consequence.loftee_ctx = Some(Box::new(ctx));
        consequence.exon_ids = Some(exon_ids);
        consequence
    }

    fn span(start: u64, end: u64) -> InputVariant {
        InputVariant::new("21".into(), start, end, b"C".to_vec(), b"T".to_vec())
    }

    fn plugin(params: &[&str]) -> NearestExonJbPlugin {
        let mut plugin = NearestExonJbPlugin::new();
        plugin
            .init(&params.iter().map(|p| p.to_string()).collect::<Vec<_>>())
            .unwrap();
        plugin
    }

    fn value(
        plugin: &NearestExonJbPlugin,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Option<String> {
        plugin
            .run(consequence, variant)
            .unwrap()
            .get(OUTPUT_FIELD)
            .cloned()
    }

    /// Three exons on the forward strand: 1000-1100, 2000-2100, 3000-3100.
    const FORWARD: [(u64, u64); 3] = [(1000, 1100), (2000, 2100), (3000, 3100)];

    /// Ensembl VEP writes `ENSE00000000002|10|start|101` for a substitution at
    /// 2010 (10 from the exon's start, 90 from its end) and `|5|end|` for one
    /// at 2095; at 2050 both junctions are 50 away: `start_end`.
    #[test]
    fn an_exonic_variant_reports_its_nearer_junction() {
        let plugin = plugin(&[]);
        let tc = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        assert_eq!(
            value(&plugin, &tc, &span(2010, 2010)).as_deref(),
            Some("ENSE00000000002|10|start|101")
        );
        assert_eq!(
            value(&plugin, &tc, &span(2095, 2095)).as_deref(),
            Some("ENSE00000000002|5|end|101")
        );
        assert_eq!(
            value(&plugin, &tc, &span(2050, 2050)).as_deref(),
            Some("ENSE00000000002|50|start_end|101")
        );
    }

    /// The distances are the variant's start to the exon's start and its end to
    /// the exon's end: a deletion of 2003 to 2098 is 3 from the start and 2
    /// from the end.
    #[test]
    fn a_spanning_variant_measures_each_end_against_its_junction() {
        let plugin = plugin(&[]);
        let tc = row(&[Consequence::InframeDeletion], None, context(1, &FORWARD));
        assert_eq!(
            value(&plugin, &tc, &span(2003, 2098)).as_deref(),
            Some("ENSE00000000002|2|end|101")
        );
    }

    /// An insertion between 1100 and 1101 is the pair `start 1101, end 1100`:
    /// it overlaps exon 1 (end 1100), 101 from the exon's start and 0 from its
    /// end.
    #[test]
    fn an_insertion_at_the_junction_is_at_distance_zero() {
        let plugin = plugin(&[]);
        let tc = row(
            &[Consequence::FrameshiftVariant],
            None,
            context(1, &FORWARD),
        );
        assert_eq!(
            value(&plugin, &tc, &span(1101, 1100)).as_deref(),
            Some("ENSE00000000001|0|end|101")
        );
    }

    /// An allele of a `--minimal` split record reads its own trimmed span, the
    /// one its row's alleles describe, for its exon set and its distances, as
    /// the same allele does when it is the record's only alternate: the `-`
    /// allele of `CAA/C/T` eight bases into an exon measures 9 from its own
    /// start, the insertion allele of `C/CA/T` there likewise, the `-` allele of
    /// `CAA/C/T` whose own span lies past the exon's last base overlaps no exon
    /// whatever its untrimmed sibling overlaps, and the `G` allele of `CA/CT/GA`
    /// on the exon's last base is at distance 0 from its end. The record's
    /// span, the one its `Location` prints, would give `|8|start|`, `|2|end|`
    /// and `|1|end|`.
    #[test]
    fn a_split_allele_reads_its_own_span() {
        let plugin = plugin(&[]);
        let tc = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        let mut deletion = span(2009, 2010);
        deletion.original_start = Some(2008);
        deletion.original_end = Some(2010);
        assert_eq!(
            value(&plugin, &tc, &deletion).as_deref(),
            Some("ENSE00000000002|9|start|101")
        );
        let mut insertion = span(2009, 2008);
        insertion.original_start = Some(2008);
        insertion.original_end = Some(2008);
        assert_eq!(
            value(&plugin, &tc, &insertion).as_deref(),
            Some("ENSE00000000002|9|start|101")
        );
        assert_eq!(
            value(&plugin, &tc, &span(2009, 2010)).as_deref(),
            Some("ENSE00000000002|9|start|101")
        );
        let intronic = row(
            &[Consequence::IntronVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        let mut straddle = span(1101, 1102);
        straddle.original_start = Some(1100);
        straddle.original_end = Some(1102);
        straddle.record_allele_string_multi = Some("CAA/C/T".into());
        assert_eq!(value(&plugin, &intronic, &straddle), None);
        let mut all_trimmed = span(1101, 1101);
        all_trimmed.original_start = Some(1100);
        all_trimmed.original_end = Some(1101);
        all_trimmed.record_allele_string_multi = Some("CA/C/CT".into());
        assert_eq!(value(&plugin, &intronic, &all_trimmed), None);
        // The alleles of `CA/CT/GA` trim to different spans: `GA` to 1100, on
        // the exon's last base; `CT` to 1101, past it.
        let mut on_exon = span(1100, 1100);
        on_exon.original_start = Some(1100);
        on_exon.original_end = Some(1101);
        on_exon.record_allele_string_multi = Some("CA/CT/GA".into());
        assert_eq!(
            value(&plugin, &tc, &on_exon).as_deref(),
            Some("ENSE00000000001|0|end|101")
        );
        let mut past_exon = span(1101, 1101);
        past_exon.original_start = Some(1100);
        past_exon.original_end = Some(1101);
        past_exon.record_allele_string_multi = Some("CA/CT/GA".into());
        assert_eq!(value(&plugin, &intronic, &past_exon), None);
        assert_eq!(value(&plugin, &intronic, &span(1101, 1102)), None);
    }

    /// A deletion across the first intron overlaps two exons: both are measured,
    /// and the exon whose junction is nearer is the one reported.
    #[test]
    fn two_overlapped_exons_report_the_nearer_junction_alone() {
        let plugin = plugin(&[]);
        let tc = row(
            &[Consequence::FrameshiftVariant],
            None,
            context(1, &FORWARD),
        );
        // 1090..2003: exon 1 is 90 from its start (1090) and 903 from its end
        // (2003 against 1100); exon 2 is 910 from its start and 97 from its
        // end: exon 1 at 90 wins.
        assert_eq!(
            value(&plugin, &tc, &span(1090, 2003)).as_deref(),
            Some("ENSE00000000001|90|start|101")
        );
        // Equal nearest distances list both exons, in transcript order.
        assert_eq!(
            value(&plugin, &tc, &span(1095, 2005)).as_deref(),
            Some("ENSE00000000001|95|start|101,ENSE00000000002|95|end|101")
        );
    }

    /// A row whose nearest junction is farther than `max_range` carries nothing;
    /// one exactly at the range is reported.
    #[test]
    fn max_range_bounds_the_reported_distance() {
        let long = context(1, &[(1000, 50000)]);
        let tc = row(&[Consequence::MissenseVariant], None, long);
        assert_eq!(value(&plugin(&[]), &tc, &span(25000, 25000)), None);
        assert_eq!(
            value(&plugin(&["max_range=30000"]), &tc, &span(25000, 25000)).as_deref(),
            Some("ENSE00000000001|24000|start|49001")
        );
        assert_eq!(
            value(&plugin(&["max_range=24000"]), &tc, &span(25000, 25000)).as_deref(),
            Some("ENSE00000000001|24000|start|49001")
        );
        assert_eq!(
            value(&plugin(&["max_range=23999"]), &tc, &span(25000, 25000)),
            None
        );
    }

    #[test]
    fn an_intronic_or_flanking_row_carries_nothing_by_default() {
        let plugin = plugin(&[]);
        let tc = row(
            &[Consequence::IntronVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        assert_eq!(value(&plugin, &tc, &span(1500, 1500)), None);
        let upstream = row(
            &[Consequence::UpstreamGeneVariant],
            None,
            context(1, &FORWARD),
        );
        assert_eq!(value(&plugin, &upstream, &span(500, 500)), None);
    }

    /// Under `intronic=1` Ensembl VEP reports the two exons flanking the intron
    /// with their own distances: at 1500 in intron 1, exon 1 (end 1100) is 400
    /// away on its end side and exon 2 (start 2000) 500 away on its start side;
    /// on the reverse strand the flanking exons are the model's exons `n` and
    /// `n + 1` all the same.
    #[test]
    fn intronic_reports_both_flanking_exons() {
        let plugin = plugin(&["intronic=1"]);
        let tc = row(
            &[Consequence::IntronVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        assert_eq!(
            value(&plugin, &tc, &span(1500, 1500)).as_deref(),
            Some("ENSE00000000001|400|end|101,ENSE00000000002|500|start|101")
        );
        let reverse = row(
            &[Consequence::IntronVariant],
            Some("2/2"),
            context(-1, &[(3000, 3100), (2000, 2100), (1000, 1100)]),
        );
        // Intron 2 of the reverse transcript lies between exon 2 (2000-2100)
        // and exon 3 (1000-1100): 1500 is 500 from exon 2's start and 400
        // from exon 3's end.
        assert_eq!(
            value(&plugin, &reverse, &span(1500, 1500)).as_deref(),
            Some("ENSE00000000002|500|start|101,ENSE00000000003|400|end|101")
        );
    }

    /// The intronic path needs an `intron_variant` term and an `INTRON` number;
    /// a variant spanning two introns reads the column's first number, and
    /// `max_range` does not apply to the flanking exons.
    #[test]
    fn intronic_reads_the_intron_number_and_ignores_max_range() {
        let plugin = plugin(&["intronic=1", "max_range=10"]);
        let tc = row(
            &[Consequence::IntronVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        assert_eq!(
            value(&plugin, &tc, &span(1500, 1500)).as_deref(),
            Some("ENSE00000000001|400|end|101,ENSE00000000002|500|start|101")
        );
        let two_introns = row(
            &[Consequence::IntronVariant],
            Some("1-2/2"),
            context(1, &FORWARD),
        );
        assert_eq!(
            value(&plugin, &two_introns, &span(1500, 1500)).as_deref(),
            Some("ENSE00000000001|400|end|101,ENSE00000000002|500|start|101")
        );
        let no_term = row(
            &[Consequence::UpstreamGeneVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        assert_eq!(value(&plugin, &no_term, &span(1500, 1500)), None);
        let no_number = row(&[Consequence::IntronVariant], None, context(1, &FORWARD));
        assert_eq!(value(&plugin, &no_number, &span(1500, 1500)), None);
    }

    /// Under `intronic=1` two overlapped exons take the flanking branch too:
    /// both are reported with their own distances (`NearestExonJB.pm:244`).
    #[test]
    fn intronic_reports_two_overlapped_exons_with_their_own_distances() {
        let plugin = plugin(&["intronic=1"]);
        let tc = row(
            &[Consequence::FrameshiftVariant],
            None,
            context(1, &FORWARD),
        );
        assert_eq!(
            value(&plugin, &tc, &span(1090, 2003)).as_deref(),
            Some("ENSE00000000001|90|start|101,ENSE00000000002|97|end|101")
        );
    }

    /// Two alleles of one record share the span 1101-1104, the first four bases
    /// of intron 1, and differ in their terms: the `splice_donor_variant` row
    /// carries no `intron_variant`, so the intronic branch reports nothing on
    /// it, while its sibling's row gets both flanking exons (exon 1 is 4 from
    /// its end against 101 from its start, exon 2 is 899 from its start against
    /// 996 from its end). Ensembl VEP caches one string per transcript and span
    /// (`NearestExonJB.pm:179`) and hands the first row's items to the second;
    /// each row here answers for its own terms, whichever runs first.
    #[test]
    fn each_row_answers_for_its_own_terms_under_intronic() {
        let plugin = plugin(&["intronic=1"]);
        let donor = row(
            &[Consequence::SpliceDonorVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        let region = row(
            &[
                Consequence::SpliceDonorRegionVariant,
                Consequence::IntronVariant,
            ],
            Some("1/2"),
            context(1, &FORWARD),
        );
        let shared = span(1101, 1104);
        let flanking = "ENSE00000000001|4|end|101,ENSE00000000002|899|start|101";
        assert_eq!(value(&plugin, &region, &shared).as_deref(), Some(flanking));
        assert_eq!(value(&plugin, &donor, &shared), None);
        assert_eq!(value(&plugin, &region, &shared).as_deref(), Some(flanking));
        assert_eq!(value(&plugin, &donor, &shared), None);
    }

    /// A transcript with an intron of at most 13 bases stretches every exon by
    /// 12 bases in the overlap test, so a variant in such an intron overlaps
    /// both of its exons.
    #[test]
    fn a_frameshift_intron_stretches_the_exons() {
        let plugin = plugin(&[]);
        let tc = row(
            &[Consequence::IntronVariant],
            Some("1/1"),
            context(1, &[(1000, 1100), (1110, 1200)]),
        );
        assert_eq!(
            value(&plugin, &tc, &span(1105, 1105)).as_deref(),
            Some("ENSE00000000001|5|end|101,ENSE00000000002|5|start|91")
        );
    }

    /// Abutting exons leave an intron whose start is one past its end, and
    /// overlapping exons one whose start is further past it (`Intron->new`);
    /// Ensembl VEP skips such an intron before its frameshift test
    /// (`BaseTranscriptVariation.pm:1014, :1032`), so neither transcript
    /// stretches its exons: a variant 5 bases past the last exon overlaps
    /// nothing, one on the exon's last base is at its end.
    #[test]
    fn an_intron_whose_start_is_past_its_end_stretches_no_exon() {
        let plugin = plugin(&[]);
        let abutting = ([(1000, 1100), (1101, 1200)], "ENSE00000000002|0|end|100");
        let overlapping = ([(1000, 1100), (1095, 1200)], "ENSE00000000002|0|end|106");
        for (exons, on_the_last_base) in [abutting, overlapping] {
            let tc = row(
                &[Consequence::DownstreamGeneVariant],
                None,
                context(1, &exons),
            );
            assert_eq!(value(&plugin, &tc, &span(1205, 1205)), None, "{exons:?}");
            assert_eq!(
                value(&plugin, &tc, &span(1200, 1200)).as_deref(),
                Some(on_the_last_base),
                "{exons:?}"
            );
        }
    }

    #[test]
    fn the_vcf_format_joins_the_parts_with_plus() {
        let mut plugin = plugin(&[]);
        plugin.set_output_format(OutputFormat::Vcf);
        let tc = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        assert_eq!(
            value(&plugin, &tc, &span(2010, 2010)).as_deref(),
            Some("ENSE00000000002+10+start+101")
        );
        assert_eq!(
            plugin.header_info()[0].1,
            "Nearest Exon Junction Boundary (coding sequence variants only). \
             Format:ExonID+distance+start/end+length"
        );
        plugin.set_output_format(OutputFormat::Tab);
        plugin.init(&["intronic=1".to_string()]).unwrap();
        assert_eq!(
            plugin.header_info()[0].1,
            "Nearest Exon Junction Boundary. Format:ExonID|distance|start/end|length"
        );
    }

    /// Under `--shift_3prime 1` a row's terms describe the shifted span while the
    /// plugin's distances are the record's (`NearestExonJB.pm:225-226`), so the
    /// exon set is read at the record's span too: a deletion of an exon's last
    /// base whose row shifted one base into the intron reports that exon at
    /// distance 0, as both engines write for the same record without the flag,
    /// and a deletion of the intron base before an exon whose row shifted onto
    /// the exon reports nothing, under `intronic=1` as well. Ensembl VEP's
    /// plugin reads the exon set its consequence tests memoised at the shifted
    /// span, the documented class `nearest_exon_jb_shifted_span`.
    #[test]
    fn a_shifted_indel_reads_the_exon_set_at_the_record_span() {
        let mut off_the_exon = row(
            &[Consequence::SpliceDonorVariant],
            None,
            context(1, &FORWARD),
        );
        off_the_exon.shift_length = 1;
        assert_eq!(
            value(&plugin(&[]), &off_the_exon, &span(1100, 1100)).as_deref(),
            Some("ENSE00000000001|0|end|101")
        );
        let mut onto_the_exon = row(
            &[
                Consequence::FrameshiftVariant,
                Consequence::SpliceRegionVariant,
            ],
            None,
            context(1, &FORWARD),
        );
        onto_the_exon.shift_length = 1;
        assert_eq!(value(&plugin(&[]), &onto_the_exon, &span(1999, 1999)), None);
        assert_eq!(
            value(&plugin(&["intronic=1"]), &onto_the_exon, &span(1999, 1999)),
            None
        );
    }

    #[test]
    fn rows_without_the_model_or_of_structural_variants_carry_nothing() {
        let plugin = plugin(&[]);
        let mut bare = TranscriptConsequence::default();
        bare.feature_type = FeatureType::Transcript;
        assert_eq!(value(&plugin, &bare, &span(2010, 2010)), None);
        let tc = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        let mut deletion = span(2010, 2500);
        deletion.is_structural = true;
        assert_eq!(value(&plugin, &tc, &deletion), None);
    }

    /// Ensembl VEP writes an exon that has no stable id with the empty id,
    /// `|<distance>|<side>|<length>` (`NearestExonJB.pm:220, :269`): on a RefSeq
    /// aligned transcript under `--all_refseq` a substitution 23 bases into a
    /// 168-base exon is `|23|start|168`. A model whose exons all lack ids
    /// carries no id list; a model with some reports the ids it has. Ensembl
    /// VEP keys the exons of a row by id, so two id-less exons share one key
    /// and it writes one item from the parts of both (`:228-237`); each exon
    /// keeps its own item here, in the default and in the intronic branch.
    #[test]
    fn an_exon_without_an_id_is_written_with_the_empty_id() {
        let default = plugin(&[]);
        let mut without_ids = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        without_ids.exon_ids = None;
        assert_eq!(
            value(&default, &without_ids, &span(2010, 2010)).as_deref(),
            Some("|10|start|101")
        );
        assert_eq!(
            value(&default, &without_ids, &span(1095, 2005)).as_deref(),
            Some("|95|start|101,|95|end|101")
        );
        let mut some_ids = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        some_ids.exon_ids = Some(Arc::from([
            Arc::from("ENSE00000000001"),
            Arc::from(""),
            Arc::from("ENSE00000000003"),
        ]));
        assert_eq!(
            value(&default, &some_ids, &span(2010, 2010)).as_deref(),
            Some("|10|start|101")
        );
        assert_eq!(
            value(&default, &some_ids, &span(1010, 1010)).as_deref(),
            Some("ENSE00000000001|10|start|101")
        );
        let mut intronic = row(
            &[Consequence::IntronVariant],
            Some("1/2"),
            context(1, &FORWARD),
        );
        intronic.exon_ids = None;
        assert_eq!(
            value(&plugin(&["intronic=1"]), &intronic, &span(1500, 1500)).as_deref(),
            Some("|400|end|101,|500|start|101")
        );
        let mut vcf = plugin(&[]);
        vcf.set_output_format(OutputFormat::Vcf);
        assert_eq!(
            value(&vcf, &without_ids, &span(2010, 2010)).as_deref(),
            Some("+10+start+101")
        );
    }

    /// Two id-less exons of one row each keep their own distance, side and
    /// length. Ensembl VEP files both under the empty id (`NearestExonJB.pm:220,
    /// :228-237`), so its item takes the length, and on a tie the side, of the
    /// exon filed last, and under `intronic=1` the distance of whichever exon
    /// its hash order reads last (`:251-259`). The model's exons:
    /// 10521560-10521694 (135 bases), 10524580-10524688 (109),
    /// 10578337-10578418 (82) and 10578517-10578605 (89). A deletion of
    /// 10578401-10578530 is 64 from the third exon's start and 75 from the
    /// fourth's end, so the third alone is reported with its own 82 bases, not
    /// 89; extended to 10578541 both are 64 away and both are reported; a
    /// substitution at 10521800 in the first intron reports the first exon at
    /// 106 with its own 135 bases and the second at 2780, never one exon's
    /// distance with the other's length.
    #[test]
    fn id_less_exons_keep_their_own_lengths_and_junctions() {
        let model = [
            (10521560, 10521694),
            (10524580, 10524688),
            (10578337, 10578418),
            (10578517, 10578605),
        ];
        let default = plugin(&[]);
        let mut exonic = row(
            &[Consequence::SpliceAcceptorVariant],
            None,
            context(1, &model),
        );
        exonic.exon_ids = None;
        assert_eq!(
            value(&default, &exonic, &span(10578401, 10578530)).as_deref(),
            Some("|64|start|82")
        );
        assert_eq!(
            value(&default, &exonic, &span(10578401, 10578541)).as_deref(),
            Some("|64|start|82,|64|end|89")
        );
        let mut intronic = row(
            &[Consequence::IntronVariant],
            Some("1/3"),
            context(1, &model),
        );
        intronic.exon_ids = None;
        assert_eq!(
            value(
                &plugin(&["intronic=1"]),
                &intronic,
                &span(10521800, 10521800)
            )
            .as_deref(),
            Some("|106|end|135,|2780|start|109")
        );
    }

    #[test]
    fn parameters_are_checked_at_init() {
        let message = |params: &[&str]| {
            NearestExonJbPlugin::new()
                .init(&params.iter().map(|p| p.to_string()).collect::<Vec<_>>())
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            message(&["50"]),
            "plugin initialization failed: NearestExonJB: parameter '50' is not key=value; \
             use max_range=<bases> or intronic=1"
        );
        assert!(message(&["gff3=genes.gff3.gz"]).contains("takes no file"));
        assert_eq!(
            message(&["parent_only=2"]),
            "plugin initialization failed: NearestExonJB: parent_only must be 0 or 1, got '2'"
        );
        assert!(message(&["range=5"]).contains("unknown parameter 'range'"));
        let plugin = plugin(&["max_range=500", "intronic=0"]);
        assert_eq!(
            (plugin.max_range, plugin.at_range, plugin.intronic),
            (500.0, Some(500), false)
        );
    }

    /// A `gff3` value Perl reads as false, `0` or the empty string, is no file
    /// (`NearestExonJB.pm:115`), and `parent_only` is read whether or not a
    /// file is (`:106-108`), so these shapes run on the cache's model as a
    /// bare `NearestExonJB` does.
    #[test]
    fn a_false_gff3_value_and_parent_only_leave_the_cache_model_in_use() {
        for (params, max_range) in [
            (&["gff3=0"][..], DEFAULT_MAX_RANGE),
            (&["gff3="][..], DEFAULT_MAX_RANGE),
            (&["parent_only=0"][..], DEFAULT_MAX_RANGE),
            (&["parent_only=1"][..], DEFAULT_MAX_RANGE),
            (&["gff3=0", "parent_only=1", "max_range=500"][..], 500),
        ] {
            let plugin = plugin(params);
            assert_eq!(
                (plugin.max_range, plugin.at_range, plugin.intronic),
                (max_range as f64, Some(max_range), false),
                "{params:?}"
            );
        }
    }

    /// `$CONFIG{intronic} == 1` (`NearestExonJB.pm:135, :190, :244`) numifies
    /// the value, so Ensembl VEP runs on every one of these shapes: `01`,
    /// `1.0`, `1e0`, `+1`, ` 1` and `1abc` are 1, and the empty string, `2`,
    /// `0.5`, `yes`, `0x1` and `10` are not.
    #[test]
    fn intronic_is_read_as_perl_compares_it_with_one() {
        for value in ["1", "01", "1.0", "1e0", "+1", " 1", "1abc"] {
            let plugin = plugin(&[&format!("intronic={value}")]);
            assert!(plugin.intronic, "intronic={value}");
            assert!(plugin.header_info()[0]
                .1
                .starts_with("Nearest Exon Junction Boundary. "));
        }
        for value in ["0", "", "2", "0.5", "yes", "0x1", "10"] {
            let plugin = plugin(&[&format!("intronic={value}")]);
            assert!(!plugin.intronic, "intronic={value}");
            assert!(plugin.header_info()[0]
                .1
                .contains("(coding sequence variants only)"));
        }
    }

    /// `$min` starts as the value's text and `$min > $startD` numifies it
    /// (`NearestExonJB.pm:186, :230-238`), so `24000.5`, `2.5e4`, `24001abc`
    /// and `inf` admit a distance of 24000. When no exon lies nearer, `$min`
    /// keeps the text and `exists $dists{$exon}{$min}` (`:268`) matches it
    /// against the distances' decimal text: an exon exactly at the range is
    /// reported under `24000` and not under `24000.0`, `024000`, `+24000` or
    /// `2.4e4`. A value with no number, `abc` or the empty string, is 0 and
    /// admits no exon, as do `-5` and `nan`.
    #[test]
    fn max_range_is_read_as_perl_numifies_it() {
        let long = context(1, &[(1000, 50000)]);
        let tc = row(&[Consequence::MissenseVariant], None, long);
        let at = |max_range: &str| {
            let plugin = plugin(&[&format!("max_range={max_range}")]);
            value(&plugin, &tc, &span(25000, 25000))
        };
        let item = Some("ENSE00000000001|24000|start|49001");
        for admitted in ["24000", "24000.5", "2.5e4", "24001abc", "inf"] {
            assert_eq!(at(admitted).as_deref(), item, "max_range={admitted}");
        }
        for refused in [
            "24000.0", "024000", "+24000", "2.4e4", "23999.9", "abc", "", "-5", "nan",
        ] {
            assert_eq!(at(refused), None, "max_range={refused}");
        }
        let plugin = plugin(&["max_range=50.0"]);
        assert_eq!((plugin.max_range, plugin.at_range), (50.0, None));
    }

    /// A junction the variant sits on is at distance 0: `max_range=0` reports
    /// it, since `$min` keeps the text `0` that the distance's key matches
    /// (`NearestExonJB.pm:268`), while `0.0` does not. Ensembl VEP's `%dists`
    /// keys the exon's length under `len` beside its distances (`:229-237`), so
    /// `max_range=len` makes Ensembl VEP report every overlapped exon with `len`
    /// as its distance; vep-rs reads `len` as 0 and reports nothing, the
    /// documented class `nearest_exon_jb_len_key_collision`.
    #[test]
    fn a_zero_range_reports_a_junction_the_variant_sits_on() {
        let tc = row(&[Consequence::MissenseVariant], None, context(1, &FORWARD));
        assert_eq!(
            value(&plugin(&["max_range=0"]), &tc, &span(2000, 2000)).as_deref(),
            Some("ENSE00000000002|0|start|101")
        );
        assert_eq!(
            value(&plugin(&["max_range=0.0"]), &tc, &span(2000, 2000)),
            None
        );
        assert_eq!(
            value(&plugin(&["max_range=len"]), &tc, &span(2000, 2000)),
            None
        );
    }
}
