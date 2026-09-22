// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! LoFTEE plugin: Loss-Of-Function Transcript Effect Estimator.
//!
//! A native Rust reimplementation of the Ensembl LoFTEE plugin
//! (github.com/konradjk/loftee), targeting parity with the upstream `grch38`
//! branch run with its default settings (splice-prediction layer OFF). LoFTEE
//! is a *consequence-filtering* plugin: for variants predicted to be
//! loss-of-function (`stop_gained`, `frameshift_variant`,
//! `splice_acceptor_variant`, `splice_donor_variant`) on a protein-coding
//! transcript, it emits a confidence call and the reasons it is not
//! high-confidence.
//!
//! **Output fields** (Perl LoFTEE field names, for concordance parity):
//! - `LoF`: `HC` (high confidence) or `LC` (low confidence)
//! - `LoF_filter`: comma-joined reasons the call is LC (e.g. `SMALL_INTRON`)
//! - `LoF_flags`: comma-joined warning flags (e.g. `SINGLE_EXON`, `NON_CAN_SPLICE`)
//! - `LoF_info`: supporting values (e.g. `INTRON_SIZE:1200`, `PERCENTILE:0.97`)
//!
//! **Filters** (LC reasons):
//! - `END_TRUNC`: PTC near the transcript end. GERP-weighted by default (needs
//!   a GERP source: `gerp_bigwig=` on GRCh38, `gerp_tabix=` on GRCh37, or
//!   `gerp=` to pick by extension); falls back to the 5%-percentile rule when
//!   `use_gerp_end_trunc=false`.
//! - `SMALL_INTRON`: splice variant in an intron below `min_intron_size` (15 bp).
//! - `GC_TO_GT_DONOR`: donor variant strengthening a GC donor toward canonical GT.
//! - `5UTR_SPLICE` / `3UTR_SPLICE`: essential-splice LoF that only affects a UTR.
//! - `ANC_ALLELE`: the LoF allele is the primate ancestral state (SNP-only;
//!   needs `human_ancestor_fa`).
//! - `INCOMPLETE_CDS`: `cds_start_NF`/`cds_end_NF` (OFF by default, like Perl).
//! - `EXON_INTRON_UNDEF`: exon/intron boundaries undefined for the variant.
//!
//! **Flags** (warnings, do not downgrade HC):
//! - `SINGLE_EXON`, `NON_CAN_SPLICE`, `NAGNAG_SITE`, `NO_EXON_NUMBER`.
//!   (`PHYLOCSF_WEAK` / `PHYLOCSF_UNLIKELY_ORF` are not emitted, see below.)
//!
//! **Not implemented:** the MaxEntScan + logistic-regression + de-novo-donor-SVM
//! splice-prediction layer and the `OS` confidence class (both off by default on the
//! upstream `grch38` branch), and the PhyloCSF flags: `conservation_file=` is parsed and ignored.
//!
//! **Usage:** `--plugin LoFTEE[,human_ancestor_fa=PATH][,gerp_bigwig=PATH]`
//! `[,gerp_tabix=PATH][,gerp=PATH][,conservation_file=PATH][,min_intron_size=15]`
//! `[,use_gerp_end_trunc=true][,gerp_end_trunc_cutoff=N][,filter_position=0.05][,check_complete_cds=false]`
//!
//! **Assembly handling.** The two upstream LoFTEE data bundles are structurally
//! different, so both are supported:
//!
//! | | GRCh37 (`master`) | GRCh38 (`grch38`) |
//! |---|---|---|
//! | GERP data | `GERP_scores.final.sorted.txt.gz` (tabix, one row per base) | `gerp_conservation_scores...bw` (bigWig) |
//! | parameter | `gerp_tabix=` | `gerp_bigwig=` |
//! | END_TRUNC cutoff | `+180` (`gerp_dist <= 180`) | `-58` |
//!
//! The cutoff differs in sign between branches, so it is chosen from
//! `--assembly` unless `gerp_end_trunc_cutoff=` is passed explicitly.
//!
//! The reference FASTA (`--fasta`) is required for the splice-motif filters
//! (`NON_CAN_SPLICE`, `GC_TO_GT_DONOR`); without it those checks are skipped.

use std::sync::Arc;

use indexmap::IndexMap;
use tracing::warn;

use crate::tabix::{TabixAnnotator, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::{Consequence, LofteeContext, TranscriptConsequence};
use vep_core::variant::{InputVariant, VariantClass};

const FIELD_LOF: &str = "LoF";
const FIELD_FILTER: &str = "LoF_filter";
const FIELD_FLAGS: &str = "LoF_flags";
const FIELD_INFO: &str = "LoF_info";

/// 0-based chromosome column in the GRCh37 GERP TSV (`#chrom pos GERP`).
const GERP_TSV_CHR_COL: usize = 0;
/// 0-based position column in the GRCh37 GERP TSV.
const GERP_TSV_POS_COL: usize = 1;
/// 0-based score column in the GRCh37 GERP TSV. Upstream master sums this
/// column's per-base values (`$gerp += $res[2]` in `gerp_dist.pl`).
const GERP_TSV_SCORE_COL: usize = 2;

/// LoFTEE default `min_intron_size` (bp): introns smaller than this fail SMALL_INTRON.
const DEFAULT_MIN_INTRON_SIZE: u64 = 15;
/// LoFTEE default `filter_position`: END_TRUNC if the PTC is in the last 5% of the CDS
/// (only used when GERP-weighted END_TRUNC is disabled).
const DEFAULT_FILTER_POSITION: f64 = 0.05;
/// LoFTEE `grch38` default GERP END_TRUNC cutoff: END_TRUNC when within 50 bp of the
/// last coding exon AND the GERP-weighted distance to the stop is <= this value.
const DEFAULT_GERP_END_TRUNC_CUTOFF: f64 = -58.0;

/// LoFTEE `master` (GRCh37) GERP END_TRUNC cutoff.
///
/// Upstream master hardcodes `push(@filters,'END_TRUNC') if ($d <= 50) &
/// ($gerp_dist <= 180);`, a positive 180, opposite in sign to the grch38
/// branch's -58. Both branches sum GERP the same way (per-base values over the
/// same exon walk), so the difference is a threshold choice, not a scale
/// difference: master's GRCh37 GERP data is predominantly positive where
/// grch38's bigWig is predominantly negative over the same regions.
///
/// Using the grch38 cutoff on GRCh37 data would make END_TRUNC essentially never
/// fire; using master's on GRCh38 data would make it always fire. The cutoff is
/// therefore selected by assembly unless the user overrides it explicitly.
const MASTER_GERP_END_TRUNC_CUTOFF: f64 = 180.0;

/// Where per-base GERP conservation scores come from.
///
/// The two upstream LoFTEE data bundles are structurally different, which is why
/// this is an enum rather than one path:
///
/// - **GRCh38** (`grch38` branch) ships `gerp_conservation_scores.homo_sapiens
///   .GRCh38.bw`, a bigWig.
/// - **GRCh37** (`master` branch) ships `GERP_scores.final.sorted.txt.gz`, a
///   tabix-indexed TSV of `#chrom pos GERP` with one row per base.
#[derive(Debug, Clone)]
pub enum GerpSource {
    /// bigWig (GRCh38 / `grch38` branch).
    BigWig(String),
    /// Tabix-indexed per-base TSV (GRCh37 / `master` branch).
    TabixTsv(String),
}

impl GerpSource {
    /// The underlying file path.
    pub fn path(&self) -> &str {
        match self {
            GerpSource::BigWig(p) | GerpSource::TabixTsv(p) => p,
        }
    }

    /// Classify a GERP path by extension.
    ///
    /// `.bw`/`.bigwig` is a bigWig; anything else is treated as the tabix TSV,
    /// since the upstream GRCh37 file is `.txt.gz`. Callers may override by
    /// choosing the variant directly.
    pub fn from_path(path: &str) -> Self {
        let lower = path.to_ascii_lowercase();
        if lower.ends_with(".bw") || lower.ends_with(".bigwig") {
            GerpSource::BigWig(path.to_string())
        } else {
            GerpSource::TabixTsv(path.to_string())
        }
    }
}

/// LoFTEE plugin (grch38-default parity, splice predictions off).
pub struct LofteePlugin {
    /// Minimum intron size; introns below this trigger SMALL_INTRON.
    min_intron_size: u64,
    /// Last-fraction-of-CDS cutoff for the percentile END_TRUNC fallback.
    filter_position: f64,
    /// Use the GERP-weighted END_TRUNC rule (grch38 default true). When false,
    /// the percentile rule is used and no GERP data is needed.
    use_gerp_end_trunc: bool,
    /// GERP-weighted-distance cutoff for END_TRUNC. Sign differs by assembly:
    /// -58 on grch38, +180 on master/GRCh37. Set from the assembly unless the
    /// user passes `gerp_end_trunc_cutoff=` explicitly.
    gerp_end_trunc_cutoff: f64,
    /// True once the user set `gerp_end_trunc_cutoff=` by hand, so the
    /// assembly-derived default does not overwrite it.
    gerp_cutoff_explicit: bool,
    /// GERP conservation data (END_TRUNC). `None` disables GERP END_TRUNC.
    /// Either a bigWig (GRCh38) or a tabix per-base TSV (GRCh37).
    gerp_source: Option<GerpSource>,
    /// Target assembly, injected before `init()`. Selects the END_TRUNC cutoff
    /// and, for a `gerp=` path with an ambiguous extension, the reader.
    assembly: Option<vep_core::assembly::Assembly>,
    /// Path to `human_ancestor.fa` (ANC_ALLELE). `None` (or "false") disables
    /// the ancestral-allele filter. Must be a plain (uncompressed) indexed FASTA
    /// (`.fa` + `.fai`); `scripts/data/setup_plugin_data.sh` documents decompressing
    /// the upstream `.gz` and re-indexing with `samtools faidx`.
    human_ancestor_fa: Option<String>,
    /// PhyloCSF SQLite DB path, parsed for argument parity with the Perl plugin and never read.
    conservation_file: Option<String>,
    /// Check `cds_start_NF`/`cds_end_NF` (INCOMPLETE_CDS). Off by default, like Perl.
    check_complete_cds: bool,
    /// Reference FASTA for splice-motif checks (set from `--fasta` at init).
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    /// Indexed ancestral-allele FASTA, opened from `human_ancestor_fa` at init.
    ancestor_fasta: Option<vep_fasta::IndexedFasta>,
}

impl LofteePlugin {
    pub fn new() -> Self {
        Self {
            min_intron_size: DEFAULT_MIN_INTRON_SIZE,
            filter_position: DEFAULT_FILTER_POSITION,
            use_gerp_end_trunc: true,
            gerp_end_trunc_cutoff: DEFAULT_GERP_END_TRUNC_CUTOFF,
            gerp_cutoff_explicit: false,
            gerp_source: None,
            assembly: None,
            human_ancestor_fa: None,
            conservation_file: None,
            check_complete_cds: false,
            reference_fasta: None,
            ancestor_fasta: None,
        }
    }
}

impl Default for LofteePlugin {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse `"n/total"` (e.g. the per-consequence `exon`/`intron` field) into
/// `(n, total)`, both 1-based. Returns `None` if absent or malformed.
fn parse_rank_of_total(s: Option<&str>) -> Option<(u32, u32)> {
    let s = s?;
    let (n, total) = s.split_once('/')?;
    Some((n.trim().parse().ok()?, total.trim().parse().ok()?))
}

/// True if the consequence list contains any of the given SO terms.
fn has_consequence(tc: &TranscriptConsequence, wanted: &[Consequence]) -> bool {
    tc.consequences.iter().any(|c| wanted.contains(c))
}

impl LofteePlugin {
    /// Fetch the two intronic boundary bases (5' donor dinucleotide and 3'
    /// acceptor dinucleotide) for the intron of the given 1-based `rank`,
    /// returned 5'->3' on the transcript strand. `None` if no FASTA or the
    /// intron is not found.
    ///
    /// LoFTEE reads `intron->seq` (already strand-oriented); this fetches the
    /// genomic boundary bases and reverse-complements on the reverse strand.
    fn intron_motif(
        &self,
        ctx: &LofteeContext,
        chr: &str,
        intron_rank: u32,
    ) -> Option<([u8; 2], [u8; 2])> {
        let fasta = self.reference_fasta.as_ref()?;
        let intron = ctx.introns.iter().find(|i| i.rank == intron_rank)?;
        let left = fasta.sequence(chr, intron.start, intron.start + 1)?;
        let right = fasta.sequence(chr, intron.end - 1, intron.end)?;
        if left.len() != 2 || right.len() != 2 {
            return None;
        }
        let left = [left[0].to_ascii_uppercase(), left[1].to_ascii_uppercase()];
        let right = [right[0].to_ascii_uppercase(), right[1].to_ascii_uppercase()];
        if ctx.strand == -1 {
            Some((revcomp2(right), revcomp2(left)))
        } else {
            Some((left, right))
        }
    }

    /// NON_CAN_SPLICE: intron does not start with GT and end with AG (on the
    /// transcript strand). Returns `None` when the motif can't be read (no FASTA).
    fn non_canonical_intron(
        &self,
        ctx: &LofteeContext,
        chr: &str,
        intron_rank: u32,
    ) -> Option<bool> {
        let (donor, acceptor) = self.intron_motif(ctx, chr, intron_rank)?;
        Some(&donor != b"GT" || &acceptor != b"AG")
    }

    /// GC_TO_GT_DONOR: a donor SNV where the intron starts with GC and the
    /// variant changes C->T, making the donor more canonical (GT). LoFTEE checks
    /// `feature_seq` ref C -> alt T against a GC intron start.
    fn gc_to_gt_donor(
        &self,
        ctx: &LofteeContext,
        variant: &InputVariant,
        chr: &str,
        intron_rank: u32,
    ) -> Option<bool> {
        let (donor, _) = self.intron_motif(ctx, chr, intron_rank)?;
        if &donor != b"GC" {
            return Some(false);
        }
        if variant.variant_class != VariantClass::Snv {
            return Some(false);
        }
        let (mut r, mut a) = (
            variant.ref_allele.first().copied()?.to_ascii_uppercase(),
            variant.alt_allele().first().copied()?.to_ascii_uppercase(),
        );
        if ctx.strand == -1 {
            r = complement(r);
            a = complement(a);
        }
        Some(r == b'C' && a == b'T')
    }

    /// Confirm the GERP source is present and openable, so a broken path fails
    /// at init instead of silently disabling END_TRUNC for the whole run.
    fn validate_gerp_source(source: &GerpSource) -> Result<(), PluginError> {
        let path = source.path();
        if !std::path::Path::new(path).exists() {
            return Err(PluginError::Init(format!(
                "LoFTEE: GERP data file not found: '{path}'"
            )));
        }
        match source {
            GerpSource::BigWig(p) => {
                bigtools::BigWigRead::open_file(p).map_err(|e| {
                    PluginError::Init(format!("LoFTEE: failed to open GERP bigWig '{p}': {e}"))
                })?;
            }
            GerpSource::TabixTsv(p) => {
                // Without its tabix index every interval query returns nothing.
                let config = TabixConfig {
                    file_path: p.into(),
                    chr_col: GERP_TSV_CHR_COL,
                    start_col: GERP_TSV_POS_COL,
                    end_col: None,
                    ref_col: None,
                    alt_col: None,
                    zero_based: false,
                };
                TabixAnnotator::open(config).map_err(|e| {
                    PluginError::Init(format!(
                        "LoFTEE: failed to open GERP tabix TSV '{p}': {e}. \
                         Expected the upstream GRCh37 GERP_scores.final.sorted.txt.gz \
                         with a sibling .tbi."
                    ))
                })?;
            }
        }
        Ok(())
    }
}

impl BuiltinPlugin for LofteePlugin {
    fn name(&self) -> &str {
        "LoFTEE"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                FIELD_LOF.to_string(),
                "Loss-of-function annotation (HC = High Confidence; LC = Low Confidence)"
                    .to_string(),
            ),
            (
                FIELD_FILTER.to_string(),
                "Reason for LoF not being HC".to_string(),
            ),
            (
                FIELD_FLAGS.to_string(),
                "Possible warning flags for LoF".to_string(),
            ),
            (
                FIELD_INFO.to_string(),
                "Info used for LoF annotation".to_string(),
            ),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        for param in params {
            let (key, value) = match param.split_once('=') {
                Some((k, v)) => (k.trim(), v.trim()),
                None => {
                    // LoFTEE's native argument syntax is `key:value`.
                    match param.split_once(':') {
                        Some((k, v)) => (k.trim(), v.trim()),
                        None => continue,
                    }
                }
            };
            match key {
                "min_intron_size" => {
                    self.min_intron_size = value.parse().map_err(|_| {
                        PluginError::Init(format!("LoFTEE: bad min_intron_size '{value}'"))
                    })?;
                }
                "filter_position" => {
                    self.filter_position = value.parse().map_err(|_| {
                        PluginError::Init(format!("LoFTEE: bad filter_position '{value}'"))
                    })?;
                }
                "gerp_end_trunc_cutoff" => {
                    self.gerp_end_trunc_cutoff = value.parse().map_err(|_| {
                        PluginError::Init(format!("LoFTEE: bad gerp_end_trunc_cutoff '{value}'"))
                    })?;
                    self.gerp_cutoff_explicit = true;
                }
                "use_gerp_end_trunc" => {
                    self.use_gerp_end_trunc = !matches!(value, "false" | "0" | "False");
                }
                "check_complete_cds" => {
                    self.check_complete_cds = matches!(value, "true" | "1" | "True");
                }
                "gerp_bigwig" => {
                    self.gerp_source = (value != "false" && !value.is_empty())
                        .then(|| GerpSource::BigWig(value.to_string()));
                }
                "gerp_tabix" | "gerp_scores" => {
                    self.gerp_source = (value != "false" && !value.is_empty())
                        .then(|| GerpSource::TabixTsv(value.to_string()));
                }
                "gerp" => {
                    self.gerp_source = (value != "false" && !value.is_empty())
                        .then(|| GerpSource::from_path(value));
                }
                "human_ancestor_fa" => {
                    self.human_ancestor_fa =
                        (value != "false" && !value.is_empty()).then(|| value.to_string());
                }
                "conservation_file" => {
                    self.conservation_file =
                        (value != "false" && !value.is_empty()).then(|| value.to_string());
                }
                "loftee_path" => { /* accepted for arg-string parity; unused natively */ }
                other => {
                    warn!("LoFTEE: ignoring unknown parameter '{other}'");
                }
            }
        }

        if !self.gerp_cutoff_explicit {
            if let Some(assembly) = self.assembly {
                self.gerp_end_trunc_cutoff = match assembly {
                    vep_core::assembly::Assembly::Grch37 => MASTER_GERP_END_TRUNC_CUTOFF,
                    vep_core::assembly::Assembly::Grch38 => DEFAULT_GERP_END_TRUNC_CUTOFF,
                };
            }
        }

        if self.use_gerp_end_trunc {
            if let Some(ref source) = self.gerp_source {
                Self::validate_gerp_source(source)?;
            }
        }

        if let Some(ref path) = self.human_ancestor_fa {
            let fasta = vep_fasta::IndexedFasta::from_path(path).map_err(|e| {
                PluginError::Init(format!(
                    "LoFTEE: failed to open human_ancestor_fa '{path}': {e}. \
                     Expected a plain (uncompressed) indexed FASTA with a sibling .fai."
                ))
            })?;
            self.ancestor_fasta = Some(fasta);
        }

        Ok(())
    }

    fn set_reference_fasta(&mut self, fasta: Option<Arc<vep_fasta::IndexedFasta>>) {
        self.reference_fasta = fasta;
    }

    fn set_assembly(&mut self, assembly: Option<vep_core::assembly::Assembly>) {
        self.assembly = assembly;
    }

    fn run(
        &self,
        tc: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut out = IndexMap::new();

        // Only protein-coding transcripts get a LoF call (`LoF.pm` `run`).
        if tc.biotype.as_deref() != Some("protein_coding") {
            return Ok(out);
        }

        // Without structural context, emit nothing rather than a wrong call.
        let Some(ctx) = tc.loftee_ctx.as_deref() else {
            return Ok(out);
        };

        let other_lof = has_consequence(
            tc,
            &[Consequence::StopGained, Consequence::FrameshiftVariant],
        );
        let vep_splice_lof = has_consequence(
            tc,
            &[
                Consequence::SpliceAcceptorVariant,
                Consequence::SpliceDonorVariant,
            ],
        );

        // With the splice-prediction layer off (grch38 default), confidence is
        // HC iff the variant is an essential-splice or stop/frameshift LoF.
        if !(other_lof || vep_splice_lof) {
            return Ok(out);
        }

        let chr = &variant.chr;
        let mut filters: Vec<String> = Vec::new();
        let mut flags: Vec<String> = Vec::new();
        let mut info: Vec<String> = Vec::new();

        let exon_no = parse_rank_of_total(tc.exon.as_deref());
        let intron_no = parse_rank_of_total(tc.intron.as_deref());

        let five_utr = is_utr_splice(ctx, variant, true);
        let three_utr = is_utr_splice(ctx, variant, false);

        if other_lof && ctx.cdna_coding_start.is_some() && ctx.cdna_coding_end.is_some() {
            if let Some(pct) = cds_percentile(tc, ctx) {
                info.push(format!("PERCENTILE:{pct:.6}"));
                if !self.use_gerp_end_trunc {
                    if pct >= 1.0 - self.filter_position {
                        filters.push("END_TRUNC".to_string());
                    }
                } else if let Some(ref gerp) = self.gerp_source {
                    if let Some((gerp_dist, bp_dist)) =
                        gerp_weighted_dist(ctx, chr, variant.start, gerp)
                    {
                        info.push(format!("GERP_DIST:{gerp_dist:.3}"));
                        info.push(format!("BP_DIST:{bp_dist}"));
                        if let Some((rank, _)) = exon_no {
                            let last_exon_len = last_exon_coding_length(ctx, rank);
                            let d = bp_dist as i64 - last_exon_len;
                            info.push(format!("DIST_FROM_LAST_EXON:{d}"));
                            info.push(format!(
                                "50_BP_RULE:{}",
                                if d <= 50 { "FAIL" } else { "PASS" }
                            ));
                            if d <= 50 && gerp_dist <= self.gerp_end_trunc_cutoff {
                                filters.push("END_TRUNC".to_string());
                            }
                        } else {
                            flags.push("NO_EXON_NUMBER".to_string());
                        }
                    }
                }
            }
        }

        if other_lof {
            if let Some((_, total_exons)) = exon_no {
                if total_exons == 1 {
                    flags.push("SINGLE_EXON".to_string());
                } else if self.check_complete_cds && (ctx.cds_start_nf || ctx.cds_end_nf) {
                    filters.push("INCOMPLETE_CDS".to_string());
                }
            } else {
                filters.push("EXON_INTRON_UNDEF".to_string());
            }
            if let (Some(cons), Some((rank, _))) = (self.conservation_file.as_deref(), exon_no) {
                if let Some((ann_orf, max_orf)) = phylocsf_lookup(cons, tc, rank) {
                    info.push(format!("ANN_ORF:{ann_orf}"));
                    info.push(format!("MAX_ORF:{max_orf}"));
                    if ann_orf < 0.0 {
                        flags.push(
                            if max_orf > 0.0 {
                                "PHYLOCSF_UNLIKELY_ORF"
                            } else {
                                "PHYLOCSF_WEAK"
                            }
                            .to_string(),
                        );
                    }
                }
            }
        }

        if let Some((intron_rank, _)) = intron_no {
            if ctx.introns.iter().all(|i| i.rank != intron_rank) {
                filters.push("EXON_INTRON_UNDEF".to_string());
            } else {
                if let Some(intron) = ctx.introns.iter().find(|i| i.rank == intron_rank) {
                    let intron_size = intron.end.saturating_sub(intron.start) + 1;
                    info.push(format!("INTRON_SIZE:{intron_size}"));
                    if intron_size < self.min_intron_size {
                        filters.push("SMALL_INTRON".to_string());
                    }
                }
                if vep_splice_lof {
                    if has_consequence(tc, &[Consequence::SpliceDonorVariant]) {
                        if let Some(true) = self.gc_to_gt_donor(ctx, variant, chr, intron_rank) {
                            filters.push("GC_TO_GT_DONOR".to_string());
                        }
                    }
                    if let Some(true) = self.non_canonical_intron(ctx, chr, intron_rank) {
                        flags.push("NON_CAN_SPLICE".to_string());
                    }
                    if five_utr {
                        filters.push("5UTR_SPLICE".to_string());
                    }
                    if three_utr {
                        filters.push("3UTR_SPLICE".to_string());
                    }
                }
                if has_consequence(tc, &[Consequence::SpliceAcceptorVariant])
                    && self.nagnag_site(variant, chr, ctx)
                {
                    flags.push("NAGNAG_SITE".to_string());
                }
            }
        }

        if self.ancestor_fasta.is_some() && self.is_ancestral_allele(variant) {
            filters.push("ANC_ALLELE".to_string());
        }

        let confidence = if filters.is_empty() { "HC" } else { "LC" };
        out.insert(FIELD_LOF.to_string(), confidence.to_string());
        if !filters.is_empty() {
            out.insert(FIELD_FILTER.to_string(), filters.join(","));
        }
        if !flags.is_empty() {
            out.insert(FIELD_FLAGS.to_string(), flags.join(","));
        }
        if !info.is_empty() {
            out.insert(FIELD_INFO.to_string(), info.join(","));
        }
        Ok(out)
    }
}

impl LofteePlugin {
    /// NAGNAG_SITE: SNP acceptor variant whose ±4 genomic context (9 bp,
    /// strand-oriented) matches `AG.AG`. LoFTEE only considers 9-bp SNP context.
    fn nagnag_site(&self, variant: &InputVariant, chr: &str, ctx: &LofteeContext) -> bool {
        if variant.variant_class != VariantClass::Snv {
            return false;
        }
        let Some(fasta) = self.reference_fasta.as_ref() else {
            return false;
        };
        let Some(seq) = fasta.sequence(chr, variant.start - 4, variant.start + 4) else {
            return false;
        };
        if seq.len() != 9 {
            return false;
        }
        let oriented = if ctx.strand == -1 { revcomp(&seq) } else { seq };
        oriented.len() == 9
            && oriented[0].eq_ignore_ascii_case(&b'A')
            && oriented[1].eq_ignore_ascii_case(&b'G')
            && oriented[3].eq_ignore_ascii_case(&b'A')
            && oriented[4].eq_ignore_ascii_case(&b'G')
    }

    /// ANC_ALLELE: SNP whose alt allele equals the primate ancestral base
    /// (from `human_ancestor.fa`). SNP-only, mirroring LoFTEE's
    /// `check_for_ancestral_allele` (which ignores indels and multi-base alleles).
    ///
    /// LoFTEE compares the *variant_feature_seq* (the alt on the transcript
    /// strand) to the ancestral base; the ancestral FASTA is forward-strand
    /// genomic, so the comparison uses the forward-strand alt allele. `variant`
    /// alleles are stored forward-strand, so no strand flip is needed here.
    fn is_ancestral_allele(&self, variant: &InputVariant) -> bool {
        if variant.variant_class != VariantClass::Snv {
            return false;
        }
        let Some(fasta) = self.ancestor_fasta.as_ref() else {
            return false;
        };
        match fasta.base(&variant.chr, variant.start) {
            Some(anc) => variant
                .alt_allele()
                .first()
                .is_some_and(|a| a.eq_ignore_ascii_case(&anc)),
            None => false,
        }
    }
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

fn revcomp2(d: [u8; 2]) -> [u8; 2] {
    [complement(d[1]), complement(d[0])]
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// CDS percentile of the variant's CDS end position (LoF.pm `get_position`,
/// fast mode): `cds_end / cds_length`. Uses the per-consequence `cds_position`
/// and the cDNA coding span from the context.
fn cds_percentile(tc: &TranscriptConsequence, ctx: &LofteeContext) -> Option<f64> {
    let cds_len = ctx.cdna_coding_end? as i64 - ctx.cdna_coding_start? as i64 + 1;
    if cds_len <= 0 {
        return None;
    }
    // cds_position may be "n" or "n-m"; take the max (the CDS end).
    let cds_pos = tc.cds_position.as_deref()?;
    let end = cds_pos
        .split(['-', '/'])
        .filter_map(|p| p.trim().parse::<i64>().ok())
        .max()?;
    Some(end as f64 / cds_len as f64)
}

/// 5UTR_SPLICE / 3UTR_SPLICE: an essential-splice variant whose genomic
/// position is outside the CDS on the 5' (`five`) or 3' side (LoF.pm
/// utr_splice.pl `check_5UTR` / `check_3UTR`).
fn is_utr_splice(ctx: &LofteeContext, variant: &InputVariant, five: bool) -> bool {
    let (Some(cds_start), Some(cds_end)) = (ctx.coding_region_start, ctx.coding_region_end) else {
        return false;
    };
    let pos = variant.start;
    if five {
        if ctx.strand == 1 {
            variant.end < cds_start
        } else {
            pos > cds_end
        }
    } else if ctx.strand == 1 {
        pos > cds_end
    } else {
        variant.end < cds_start
    }
}

/// Coding length of the last CDS exon at/after the variant exon (LoF.pm
/// `get_last_exon_coding_length`). Returns a bp length, or a large sentinel if
/// it cannot be located (matches LoFTEE's `-1000`, which makes `d` large so the
/// 50-bp rule does not fire).
fn last_exon_coding_length(ctx: &LofteeContext, _variant_exon_rank: u32) -> i64 {
    let stop_codon_pos = if ctx.strand == 1 {
        ctx.coding_region_end
    } else {
        ctx.coding_region_start
    };
    let Some(stop) = stop_codon_pos else {
        return -1000;
    };
    for exon in &ctx.exons {
        if ctx.strand == 1 {
            if exon.start <= stop && exon.end >= stop {
                return (stop - exon.start) as i64;
            }
        } else if exon.start <= stop && exon.end >= stop {
            return (exon.end - stop) as i64;
        }
    }
    -1000
}

/// GERP-weighted distance from the variant to the stop codon, summed over the
/// coding exon intervals between the variant and the stop codon. Port of LoFTEE
/// `gerp_dist.pl::get_gerp_weighted_dist` (grch38 / bigWig variant).
///
/// Returns `(weighted_dist, bp_dist)` where `weighted_dist` sums the per-base
/// GERP scores over each contributing interval and `bp_dist` is the summed
/// `abs(end - start)` of those intervals.
///
/// Intended divergence from Perl LoFTEE: Perl reads the bigWig's lossy
/// zoom-level summary bins (`length * binMean(score)` over `Bio::DB::BigWig`
/// summary features), whereas this sums the exact per-base values from
/// `get_interval`. On borderline loci the two differ enough to flip the
/// `END_TRUNC` decision around `gerp_end_trunc_cutoff`.
///
/// The bigWig is opened per call (END_TRUNC only fires for stop/frameshift LoF,
/// a small fraction of variants), mirroring the tabix engine's fresh-handle
/// pattern; `BigWigRead::get_interval` requires `&mut`, which a shared `&self`
/// plugin cannot hold.
fn gerp_weighted_dist(
    ctx: &LofteeContext,
    chr: &str,
    pos: u64,
    source: &GerpSource,
) -> Option<(f64, u64)> {
    let (stop_codon_pos, _start_codon_pos) = if ctx.strand == 1 {
        (ctx.coding_region_end?, ctx.coding_region_start?)
    } else {
        (ctx.coding_region_start?, ctx.coding_region_end?)
    };

    let mut summer = GerpSummer::open(source)?;

    let mut weighted_dist = 0.0_f64;
    let mut bp_dist = 0_u64;

    for exon in &ctx.exons {
        if ctx.strand == -1 {
            if pos < exon.start {
                continue;
            }
        } else if pos > exon.end {
            continue;
        }

        let last_exon = if ctx.strand == 1 {
            exon.start < stop_codon_pos && exon.end >= stop_codon_pos
        } else {
            exon.end > stop_codon_pos && exon.start <= stop_codon_pos
        };
        let in_affected_exon = pos >= exon.start && pos <= exon.end;

        let (start, end) = if last_exon {
            let start = if in_affected_exon {
                pos
            } else if ctx.strand == 1 {
                exon.start
            } else {
                exon.end
            };
            (start, stop_codon_pos)
        } else if in_affected_exon {
            let end = if ctx.strand == 1 {
                exon.end
            } else {
                exon.start
            };
            (pos, end)
        } else {
            (exon.start, exon.end)
        };

        let (lo, hi) = (start.min(end), start.max(end));
        weighted_dist += summer.sum_interval(chr, lo, hi);
        bp_dist += hi.abs_diff(lo);
    }

    Some((weighted_dist, bp_dist))
}

/// Per-format GERP interval summer.
///
/// Both upstream branches walk the same exon intervals; they differ only in how
/// they total one interval:
///
/// - **bigWig** (grch38): `sum(interval_len * value)` over the stored intervals
///   overlapping the range, equal to `length * binMean` over bigWig summary
///   features. Coordinates are 0-based, half-open `[start, end)`.
/// - **tabix TSV** (master/GRCh37): the file holds one row per base, and
///   `gerp_dist.pl` does `$gerp += $res[2]`, a plain sum of per-base values with
///   no length multiplier. Coordinates are 1-based inclusive.
///
/// Handles are opened once per `gerp_weighted_dist` call rather than per interval
/// (END_TRUNC only fires for stop/frameshift LoF, a small fraction of variants).
/// `BigWigRead::get_interval` needs `&mut`, which a shared `&self` plugin cannot
/// hold, so the reader lives here rather than on the plugin.
enum GerpSummer {
    BigWig(bigtools::BigWigRead<bigtools::utils::reopen::ReopenableFile>),
    TabixTsv(TabixAnnotator),
}

impl GerpSummer {
    fn open(source: &GerpSource) -> Option<Self> {
        match source {
            GerpSource::BigWig(path) => Some(GerpSummer::BigWig(
                bigtools::BigWigRead::open_file(path).ok()?,
            )),
            GerpSource::TabixTsv(path) => {
                let config = TabixConfig {
                    file_path: path.into(),
                    chr_col: GERP_TSV_CHR_COL,
                    start_col: GERP_TSV_POS_COL,
                    end_col: None,
                    ref_col: None,
                    alt_col: None,
                    zero_based: false,
                };
                Some(GerpSummer::TabixTsv(TabixAnnotator::open(config).ok()?))
            }
        }
    }

    /// Total GERP over the 1-based inclusive range `[lo, hi]`.
    fn sum_interval(&mut self, chr: &str, lo: u64, hi: u64) -> f64 {
        match self {
            GerpSummer::BigWig(reader) => {
                // bigWig is 0-based half-open, so shift the start and use `hi`
                // directly as the exclusive end.
                let q_start = lo.saturating_sub(1) as u32;
                let q_end = hi as u32;
                let mut total = 0.0_f64;
                if let Ok(iter) = reader.get_interval(chr, q_start, q_end) {
                    for value in iter.flatten() {
                        let clipped_start = value.start.max(q_start);
                        let clipped_end = value.end.min(q_end);
                        if clipped_end > clipped_start {
                            total += (clipped_end - clipped_start) as f64 * value.value as f64;
                        }
                    }
                }
                total
            }
            GerpSummer::TabixTsv(annotator) => {
                let mut total = 0.0_f64;
                if let Ok(records) = annotator.query(chr, lo, hi) {
                    for record in records {
                        if let Some(field) = record.columns.get(GERP_TSV_SCORE_COL) {
                            if let Ok(v) = field.trim().parse::<f64>() {
                                total += v;
                            }
                        }
                    }
                }
                total
            }
        }
    }
}

/// PhyloCSF ORF scores for the exon (LoF.pm `check_for_conservation`).
///
/// **Deliberately not implemented.** The PhyloCSF output is `PHYLOCSF_WEAK` /
/// `PHYLOCSF_UNLIKELY_ORF`, which are *flags*, not *filters*: LoFTEE adds them to
/// `LoF_flags` and they never change the HC/LC confidence call (in `LoF.pm`
/// run(), `@flags` does not feed the LC downgrade, only `@filters` does).
/// Implementing them requires a SQLite reader for the `phylocsf_summary` table,
/// a C-linked dependency, for zero HC/LC concordance benefit. The
/// `conservation_file` parameter is accepted for arg-string parity with the Perl
/// plugin and is a no-op. Returning `None` emits no PhyloCSF flag, which is also
/// what a real lookup that found no row would do.
fn phylocsf_lookup(
    _conservation_file: &str,
    _tc: &TranscriptConsequence,
    _exon_rank: u32,
) -> Option<(f64, f64)> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::ConsequenceList;

    #[test]
    fn name_and_header_fields() {
        let p = LofteePlugin::new();
        assert_eq!(p.name(), "LoFTEE");
        let fields: Vec<_> = p.header_info().into_iter().map(|(k, _)| k).collect();
        assert_eq!(fields, vec!["LoF", "LoF_filter", "LoF_flags", "LoF_info"]);
    }

    #[test]
    fn parse_rank_of_total_works() {
        assert_eq!(parse_rank_of_total(Some("3/5")), Some((3, 5)));
        assert_eq!(parse_rank_of_total(Some("1/1")), Some((1, 1)));
        assert_eq!(parse_rank_of_total(None), None);
        assert_eq!(parse_rank_of_total(Some("bad")), None);
    }

    #[test]
    fn init_parses_params() {
        // No data-file params: `init` eagerly opens `human_ancestor_fa`, and
        // `use_gerp_end_trunc=false` skips GERP-source validation.
        let mut p = LofteePlugin::new();
        p.init(&[
            "min_intron_size=20".to_string(),
            "use_gerp_end_trunc=false".to_string(),
            "gerp_end_trunc_cutoff=-60".to_string(),
            "gerp_bigwig=/path/to/loftee/grch38/gerp.bw".to_string(),
        ])
        .unwrap();
        assert_eq!(p.min_intron_size, 20);
        assert!(!p.use_gerp_end_trunc);
        assert_eq!(p.gerp_end_trunc_cutoff, -60.0);
        assert_eq!(
            p.gerp_source.as_ref().map(|s| s.path()),
            Some("/path/to/loftee/grch38/gerp.bw")
        );
    }

    /// The GRCh37 bundle ships GERP as a tabix per-base TSV, not a bigWig, so
    /// `gerp=` must classify by extension rather than assuming bigWig.
    #[test]
    fn gerp_source_classifies_by_extension() {
        assert!(matches!(
            GerpSource::from_path("/d/gerp_conservation_scores.homo_sapiens.GRCh38.bw"),
            GerpSource::BigWig(_)
        ));
        assert!(matches!(
            GerpSource::from_path("/d/GERP.BW"),
            GerpSource::BigWig(_)
        ));
        assert!(matches!(
            GerpSource::from_path("/d/GERP_scores.final.sorted.txt.gz"),
            GerpSource::TabixTsv(_)
        ));
    }

    /// `gerp_tabix=` selects the per-base TSV reader, the format the GRCh37
    /// bundle ships.
    #[test]
    fn accepts_a_tabix_gerp_source_for_grch37() {
        let mut p = LofteePlugin::new();
        p.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        p.init(&[
            "use_gerp_end_trunc=false".to_string(),
            "gerp_tabix=/path/to/loftee/grch37/GERP_scores.final.sorted.txt.gz".to_string(),
        ])
        .unwrap();
        assert!(matches!(p.gerp_source, Some(GerpSource::TabixTsv(_))));
    }

    /// The END_TRUNC cutoff differs in sign between the branches (+180 master vs
    /// -58 grch38). Applying the grch38 value to GRCh37 data would make END_TRUNC
    /// effectively never fire, so it is derived from the assembly.
    #[test]
    fn end_trunc_cutoff_follows_the_assembly() {
        let mut g37 = LofteePlugin::new();
        g37.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        g37.init(&[]).unwrap();
        assert_eq!(g37.gerp_end_trunc_cutoff, MASTER_GERP_END_TRUNC_CUTOFF);

        let mut g38 = LofteePlugin::new();
        g38.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        g38.init(&[]).unwrap();
        assert_eq!(g38.gerp_end_trunc_cutoff, DEFAULT_GERP_END_TRUNC_CUTOFF);

        let mut none = LofteePlugin::new();
        none.init(&[]).unwrap();
        assert_eq!(none.gerp_end_trunc_cutoff, DEFAULT_GERP_END_TRUNC_CUTOFF);
    }

    /// An explicit cutoff must survive assembly derivation.
    #[test]
    fn explicit_cutoff_overrides_the_assembly_default() {
        let mut p = LofteePlugin::new();
        p.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        p.init(&["gerp_end_trunc_cutoff=42.5".to_string()]).unwrap();
        assert_eq!(p.gerp_end_trunc_cutoff, 42.5);
    }

    /// A missing or unreadable GERP file must fail at init: a per-call open that
    /// returns `None` would skip END_TRUNC silently, making a run with a broken
    /// GERP path identical to a working one.
    #[test]
    fn missing_gerp_source_is_a_hard_error() {
        let mut p = LofteePlugin::new();
        let err = p
            .init(&["gerp_bigwig=/nonexistent/gerp.bw".to_string()])
            .expect_err("a missing GERP file must fail at init");
        assert!(
            err.to_string().contains("GERP data file not found"),
            "unexpected error: {err}"
        );

        let mut p2 = LofteePlugin::new();
        let err2 = p2
            .init(&["gerp_tabix=/nonexistent/GERP_scores.txt.gz".to_string()])
            .expect_err("a missing tabix GERP file must fail at init");
        assert!(err2.to_string().contains("GERP data file not found"));
    }

    /// A tabix GERP file without its `.tbi` cannot answer interval queries, so it
    /// must be rejected at init rather than silently returning zero everywhere.
    #[test]
    fn unindexed_tabix_gerp_source_is_a_hard_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("GERP_scores.txt.gz");
        std::fs::write(&path, b"not a real bgzf file").unwrap();

        let mut p = LofteePlugin::new();
        let err = p
            .init(&[format!("gerp_tabix={}", path.display())])
            .expect_err("an unindexed GERP TSV must fail at init");
        assert!(
            err.to_string().contains("failed to open GERP tabix TSV"),
            "unexpected error: {err}"
        );
    }

    /// The per-base TSV is summed verbatim (`$gerp += $res[2]` upstream), with no
    /// length multiplier, unlike the bigWig path's `length * value`.
    #[test]
    fn tabix_gerp_sums_per_base_values() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        // Real GRCh37 GERP layout: #chrom pos GERP, one row per base.
        let lines = [
            "#chrom\tpos\tGERP",
            "21\t100\t1.5",
            "21\t101\t2.0",
            "21\t102\t-0.5",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut summer =
            GerpSummer::open(&GerpSource::TabixTsv(path.to_string_lossy().to_string()))
                .expect("open tabix GERP");

        // 1.5 + 2.0 + (-0.5) = 3.0, a plain sum with no length weighting.
        let total = summer.sum_interval("21", 100, 102);
        assert!(
            (total - 3.0).abs() < 1e-9,
            "expected per-base sum 3.0, got {total}"
        );

        let partial = summer.sum_interval("21", 100, 101);
        assert!(
            (partial - 3.5).abs() < 1e-9,
            "expected 3.5 for the first two bases, got {partial}"
        );
    }

    #[test]
    fn init_missing_ancestor_fa_errors() {
        let mut p = LofteePlugin::new();
        let err = p.init(&["human_ancestor_fa=/nonexistent/human_ancestor.fa".to_string()]);
        assert!(err.is_err());
    }

    #[test]
    fn init_accepts_colon_syntax() {
        let mut p = LofteePlugin::new();
        p.init(&["min_intron_size:25".to_string()]).unwrap();
        assert_eq!(p.min_intron_size, 25);
    }

    #[test]
    fn complement_and_revcomp() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(revcomp2(*b"GT"), [b'A', b'C']); // GT donor -> AC on opposite strand
        assert_eq!(revcomp(b"GT"), b"AC");
    }

    fn coding_ctx(strand: i8) -> LofteeContext {
        LofteeContext {
            exons: Vec::new(),
            introns: Vec::new(),
            strand,
            coding_region_start: Some(1000),
            coding_region_end: Some(2000),
            translation_start: Some(1000),
            translation_end: Some(2000),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(900),
            cds_start_nf: false,
            cds_end_nf: false,
        }
    }

    fn lof_variant(
        consequences: ConsequenceList,
        ctx: LofteeContext,
    ) -> (TranscriptConsequence, InputVariant) {
        let tc = TranscriptConsequence {
            biotype: Some("protein_coding".into()),
            consequences,
            exon: Some("2/5".to_string()),
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        };
        let variant = InputVariant::new("21".into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        (tc, variant)
    }

    #[test]
    fn non_coding_biotype_returns_empty() {
        let p = LofteePlugin::new();
        let tc = TranscriptConsequence {
            biotype: Some("lncRNA".into()),
            consequences: smallvec![Consequence::StopGained],
            exon: Some("2/5".to_string()),
            loftee_ctx: Some(Box::new(coding_ctx(1))),
            ..TranscriptConsequence::default()
        };
        let v = InputVariant::new("21".into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        assert!(p.run(&tc, &v).unwrap().is_empty());
    }

    #[test]
    fn non_lof_consequence_returns_empty() {
        let p = LofteePlugin::new();
        let (tc, v) = lof_variant(smallvec![Consequence::MissenseVariant], coding_ctx(1));
        assert!(p.run(&tc, &v).unwrap().is_empty());
    }

    #[test]
    fn stop_gained_single_exon_is_hc_with_single_exon_flag() {
        let p = LofteePlugin::new();
        let tc = TranscriptConsequence {
            biotype: Some("protein_coding".into()),
            consequences: smallvec![Consequence::StopGained],
            exon: Some("1/1".to_string()),
            cds_position: Some("150".to_string()),
            loftee_ctx: Some(Box::new(coding_ctx(1))),
            ..TranscriptConsequence::default()
        };
        let v = InputVariant::new("21".into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        let out = p.run(&tc, &v).unwrap();
        assert_eq!(out.get("LoF").map(String::as_str), Some("HC"));
        assert_eq!(
            out.get("LoF_flags").map(String::as_str),
            Some("SINGLE_EXON")
        );
        assert!(out.get("LoF_filter").is_none());
    }

    #[test]
    fn small_intron_splice_is_lc() {
        let p = LofteePlugin::new();
        let mut ctx = coding_ctx(1);
        ctx.introns = vec![vep_core::transcript::Intron {
            start: 1490,
            end: 1495, // 6 bp intron < 15
            rank: 2,
        }];
        let tc = TranscriptConsequence {
            biotype: Some("protein_coding".into()),
            consequences: smallvec![Consequence::SpliceDonorVariant],
            intron: Some("2/4".to_string()),
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        };
        let v = InputVariant::new("21".into(), 1491, 1491, b"C".to_vec(), b"T".to_vec());
        let out = p.run(&tc, &v).unwrap();
        assert_eq!(out.get("LoF").map(String::as_str), Some("LC"));
        assert_eq!(
            out.get("LoF_filter").map(String::as_str),
            Some("SMALL_INTRON")
        );
        assert_eq!(
            out.get("LoF_info").map(String::as_str),
            Some("INTRON_SIZE:6")
        );
    }
}
