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
//! - `END_TRUNC`: PTC near the transcript end. GERP-weighted by default, from
//!   `gerp_bigwig=` on GRCh38, `gerp_tabix=` (or `gerp_file=`) on GRCh37, `gerp=`
//!   to pick by extension, or, when none is named, the source Ensembl VEP's
//!   plugin reads by default (the bigWig under `loftee_path` on GRCh38, master's
//!   fixed path on GRCh37); on GRCh38 it falls back to the 5%-percentile rule
//!   when `use_gerp_end_trunc=false`, and on GRCh37 the GERP rule is the only one.
//! - `SMALL_INTRON`: splice variant in an intron below `min_intron_size` (15 bp).
//! - `GC_TO_GT_DONOR`: donor variant strengthening a GC donor toward canonical GT.
//! - `5UTR_SPLICE` / `3UTR_SPLICE`: essential-splice LoF that only affects a UTR.
//! - `ANC_ALLELE`: the LoF allele is the primate ancestral state (SNP-only;
//!   needs `human_ancestor_fa`).
//! - `INCOMPLETE_CDS`: `cds_start_NF`/`cds_end_NF` (OFF by default, like Perl).
//! - `EXON_INTRON_UNDEF`: the intron number names an intron the transcript model
//!   lacks, a shape a cache built from the model never produces.
//!
//! **Flags** (warnings, do not downgrade HC):
//! - `SINGLE_EXON`, `NON_CAN_SPLICE`, `NAGNAG_SITE`, `NO_EXON_NUMBER`.
//!   (`PHYLOCSF_WEAK` / `PHYLOCSF_UNLIKELY_ORF` are not computed, see below.)
//!
//! **Not implemented:** the MaxEntScan + logistic-regression + de-novo-donor-SVM
//! splice-prediction layer and the `OS` confidence class (both off by default on the
//! upstream `grch38` branch), and the PhyloCSF items.
//! `run_splice_predictions=` is accepted so the key does not trip the unknown-parameter
//! warning when an Ensembl VEP `LoF` argument string is reused: `0` or an empty value
//! is accepted silently, and any other value asks for the layer, which is refused at
//! plugin initialization. On GRCh37, where `master` runs the layer unless the key
//! holds `0` (its `LoF.pm` lines 98 and 186), a run whose key is absent or empty is
//! refused the same way, so a GRCh37 run passes `run_splice_predictions:0`.
//! Under `conservation_file=` (off by default) `LoF.pm` reads a PhyloCSF database
//! for each stop-gained or frameshift row with an exon number and writes the
//! `LoF_info` items `ANN_ORF` and `MAX_ORF`, or `PHYLOCSF_TOO_SHORT` for an exon
//! the database lacks, and the `LoF_flags` `PHYLOCSF_WEAK` or
//! `PHYLOCSF_UNLIKELY_ORF` (lines 274-285); none is a filter, so none changes the
//! HC or LC call. No database is read here: `false` in any letter case or an
//! empty value leaves the lookup off, and any other value asks for those items
//! and is refused at plugin initialization.
//!
//! **Usage:** `--plugin LoFTEE[,loftee_path=DIR/][,human_ancestor_fa=PATH][,gerp_bigwig=PATH]`
//! `[,gerp_tabix=PATH][,gerp=PATH][,conservation_file=false][,min_intron_size=15]`
//! `[,use_gerp_end_trunc=true][,gerp_end_trunc_cutoff=N][,filter_position=0.05][,check_complete_cds=false]`
//! `[,run_splice_predictions=0]`; `LoF`, the Ensembl VEP module's name, loads the same plugin,
//! and a word holding a colon is LoFTEE's own `key:value`, read whatever `=` it
//! holds when it splits into exactly two fields. `min_intron_size`,
//! `filter_position` and `gerp_end_trunc_cutoff` are read as the numeric
//! comparisons of `LoF.pm` (lines 297, 241 and 254) read their text: the leading
//! number, 0 when none leads it.
//!
//! **Assembly handling.** The two upstream LoFTEE data bundles are structurally
//! different, so both are supported:
//!
//! | | GRCh37 (`master`) | GRCh38 (`grch38`) |
//! |---|---|---|
//! | GERP data | `GERP_scores.final.sorted.txt.gz` (tabix, one row per base) | `gerp_conservation_scores...bw` (bigWig) |
//! | parameter | `gerp_tabix=` | `gerp_bigwig=` |
//! | END_TRUNC cutoff | `+180`, a literal (`gerp_dist <= 180`) | `-58` unless `gerp_end_trunc_cutoff=` |
//! | `LoF_info` under the GERP rule | `GERP_DIST`, `BP_DIST`, `PERCENTILE` | `PERCENTILE`, `GERP_DIST`, `BP_DIST` |
//! | keys read on that branch alone | `skip_lof_info`, `tabix_path` | `use_gerp_end_trunc`, `gerp_end_trunc_cutoff`, `filter_position` |
//! | `run_splice_predictions` default | `1`: a run whose key is not `0` is refused | `0` |
//!
//! The branch is chosen from `--assembly`. On GRCh37 the three keys master
//! never reads are logged and change nothing, so the GERP rule and its cutoff
//! of 180 always apply there, a `skip_lof_info` true by Perl's string truth
//! drops `LoF_info`, and `tabix_path` is accepted and unused, the GERP TSV being
//! read in-process; on GRCh38 the last two are unknown parameters.
//!
//! The reference FASTA (`--fasta`) is required for the splice-motif filters
//! (`NON_CAN_SPLICE`, `GC_TO_GT_DONOR`); without it those checks are skipped.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{BBIReadError, BigWigRead};
use indexmap::IndexMap;
use tracing::warn;

use crate::bigwig::open_bigwig;
use crate::params::{perl_numeric, warn_unread_parameters};
use crate::perl_format::perl_number;
use crate::tabix::{TabixAnnotator, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::{Consequence, LofteeContext, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
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

/// LoFTEE default `min_intron_size` (bp, `LoF.pm` line 74): introns smaller than
/// this fail SMALL_INTRON.
const DEFAULT_MIN_INTRON_SIZE: f64 = 15.0;
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
/// therefore selected by assembly; on GRCh37 no parameter moves it, since
/// master compares with the literal.
const MASTER_GERP_END_TRUNC_CUTOFF: f64 = 180.0;

/// LoFTEE's `loftee_path` default (`LoF.pm` line 83), spelled with its trailing
/// slash because the plugin appends file names to the value as given.
const DEFAULT_LOFTEE_PATH: &str = "/vep/loftee/";
/// The GERP bigWig the `grch38` branch reads when no `gerp_bigwig` is given
/// (`LoF.pm` line 130): `loftee_path` followed by this name, with no separator
/// inserted between them.
const GRCH38_DEFAULT_GERP_BIGWIG: &str = "gerp_conservation_scores.homo_sapiens.GRCh38.bw";
/// The GERP tabix TSV the `master` branch reads when no `gerp_file` is given
/// (master `LoF.pm` line 132): a fixed path, not derived from `loftee_path`.
const GRCH37_DEFAULT_GERP_FILE: &str = "/vep/loftee/GERP_scores.final.sorted.txt.gz";
/// The note a derived GRCh38 GERP source's open failure carries after the path,
/// so the refusal says where the path came from and how to run without it.
const DERIVED_GERP_BIGWIG_NOTE: &str = "Ensembl VEP's LoF plugin derives this GERP source from \
     loftee_path when no gerp_bigwig is given; pass the GERP file, or use_gerp_end_trunc=0 for \
     the percentile rule";
/// The same note for the GRCh37 default, a fixed path `loftee_path` does not
/// move; master has no rule that runs without a GERP source.
const DERIVED_GERP_TABIX_NOTE: &str = "Ensembl VEP's LoF plugin reads this GERP source when no \
     gerp_file is given, and judges END_TRUNC by the GERP rule alone on GRCh37; pass the GERP file";

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
    /// The file path as the parameter gave it.
    pub fn path(&self) -> &str {
        match self {
            GerpSource::BigWig(p) | GerpSource::TabixTsv(p) => p,
        }
    }

    /// The file path that is opened: the bigWig path as given (the `grch38`
    /// `gerp_dist.pl` line 72 opens it by name), the tabix path without the
    /// blanks around it, since the `master` `LoF.pm` runs `tabix` on
    /// `gerp_file` as a command Perl splits into words itself (its lines 138
    /// and 141, run by `gerp_dist.pl` lines 92 and 119; the trailing `2>&1` is
    /// the one redirection Perl performs without a shell).
    pub fn opened_path(&self) -> &str {
        match self {
            GerpSource::BigWig(p) => p,
            GerpSource::TabixTsv(p) => perl_unpadded(p),
        }
    }

    /// Classify a GERP path by extension.
    ///
    /// `.bw`/`.bigwig` is a bigWig; anything else is treated as the tabix TSV,
    /// since the upstream GRCh37 file is `.txt.gz`. Callers may override by
    /// choosing the variant directly. The extension is read inside the blanks
    /// around the path, which both variants keep as given.
    pub fn from_path(path: &str) -> Self {
        let lower = perl_unpadded(path).to_ascii_lowercase();
        if lower.ends_with(".bw") || lower.ends_with(".bigwig") {
            GerpSource::BigWig(path.to_string())
        } else {
            GerpSource::TabixTsv(path.to_string())
        }
    }
}

/// LoFTEE plugin (grch38-default parity, splice predictions off).
pub struct LofteePlugin {
    /// Minimum intron size; introns below this trigger SMALL_INTRON. Perl's
    /// number of the value, as `$intron_size < $self->{min_intron_size}`
    /// (`LoF.pm` line 297) reads it: a fraction compares as given and a
    /// negative value flags no intron.
    min_intron_size: f64,
    /// Last-fraction-of-CDS cutoff for the percentile END_TRUNC fallback. Perl's
    /// number of the value, as `$lof_percentile >= 1-$self->{filter_position}`
    /// (`LoF.pm` line 241) reads it. The default on GRCh37, where master has no
    /// percentile rule and the key is not read.
    filter_position: f64,
    /// Use the GERP-weighted END_TRUNC rule (grch38 default true). When false,
    /// the percentile rule is used and no GERP data is needed. Always true on
    /// GRCh37, where master never tests the key.
    use_gerp_end_trunc: bool,
    /// GERP-weighted-distance cutoff for END_TRUNC. Sign differs by assembly:
    /// -58 on grch38, +180 on master/GRCh37. Set from the assembly; on GRCh38
    /// a `gerp_end_trunc_cutoff=` replaces it with Perl's number of its text, as
    /// `$gerp_dist <= $self->{gerp_end_trunc_cutoff}` (`LoF.pm` line 254) reads
    /// it, and on GRCh37 the key is not read, master comparing with the literal
    /// 180 (its line 268).
    gerp_end_trunc_cutoff: f64,
    /// True once `gerp_end_trunc_cutoff=` set the cutoff, so the
    /// assembly-derived default does not overwrite it.
    gerp_cutoff_explicit: bool,
    /// GERP conservation data (END_TRUNC), a bigWig (GRCh38) or a tabix per-base
    /// TSV (GRCh37). `init` fills it from the parameters or, with the GERP rule
    /// on and none named, from [`Self::default_gerp_source`]; `None` after
    /// `init` means the percentile rule.
    gerp_source: Option<GerpSource>,
    /// `gerp_source` opened at `init`, read on every row the GERP rule judges.
    gerp_reader: Option<GerpReader>,
    /// Target assembly, injected before `init()`. It picks the upstream branch
    /// followed where `master` (GRCh37) and `grch38` differ: the keys read
    /// (`use_gerp_end_trunc`, `gerp_end_trunc_cutoff` and `filter_position` on
    /// `grch38` alone, `skip_lof_info` and `tabix_path` on `master` alone), the
    /// END_TRUNC cutoff, the default GERP source and the GERP key that wins,
    /// the place of `PERCENTILE` among the GERP items of `LoF_info`, and the
    /// default of `run_splice_predictions`, which on `master` asks for the
    /// layer vep-rs refuses.
    assembly: Option<vep_core::assembly::Assembly>,
    /// Path to the ancestral-allele FASTA (ANC_ALLELE). `None` (the word `false`
    /// in any letter case, or no parameter) disables the ancestral-allele
    /// filter. An indexed FASTA, plain (`.fa` + `.fai`) or BGZF (`.fa.gz` +
    /// `.fai` + `.gzi`, the form the upstream bundle ships), opened by
    /// `vep_fasta::IndexedFasta::from_path` as `--fasta` is.
    human_ancestor_fa: Option<String>,
    /// `loftee_path` as given: the directory the GRCh38 default GERP bigWig is
    /// derived from when no GERP source is named.
    loftee_path: Option<String>,
    /// Check `cds_start_NF`/`cds_end_NF` (INCOMPLETE_CDS). Off by default, like Perl.
    check_complete_cds: bool,
    /// Write no `LoF_info`. Set on GRCh37 alone, from `skip_lof_info` by Perl's
    /// string truth: master's `LoF.pm` writes the field while the key is false
    /// (lines 243 and 339), and the grch38 branch has no such key.
    skip_lof_info: bool,
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
            gerp_reader: None,
            assembly: None,
            human_ancestor_fa: None,
            loftee_path: None,
            check_complete_cds: false,
            skip_lof_info: false,
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

/// Split one `--plugin LoF,...` word into its key and value. A word holding a
/// colon is LoFTEE's `key:value`, which `LoF.pm` (lines 67-72) reads only when
/// `split /:/` yields exactly two fields, any `=` left in the key or the value.
/// Trailing colons are trimmed first, as Perl's `split` drops the empty fields
/// they produce (`key:value:` is `key:value`); a word then holding a second
/// colon (`key:value:x`, `key::value`, the `://` of a URL) or none (`key:`,
/// `key=value:`) is dropped, as it is there. A word without a colon is vep-rs's
/// `key=value`, read at its first `=`. A word with neither is dropped, as
/// `LoF.pm` drops it. Both fields keep every blank, as `split` leaves them:
/// ` key` is another key, which `LoF.pm` stores and never reads, and ` 0` is
/// another value, a true one at its truth tests.
fn split_parameter(param: &str) -> Option<(&str, &str)> {
    if param.contains(':') {
        let (key, value) = param.trim_end_matches(':').split_once(':')?;
        if value.contains(':') {
            return None;
        }
        Some((key, value))
    } else {
        param.split_once('=')
    }
}

/// The word that switches a LoFTEE data-file parameter off: `false` in any
/// letter case, the reading of `LoF.pm`'s `lc(...) ne 'false'` tests (lines 274
/// and 312), or an empty value, which only a `key=` word gives: `split_parameter`
/// drops a bare `key:`, and a blank after the colon is a value.
fn is_off_word(value: &str) -> bool {
    value.is_empty() || value.eq_ignore_ascii_case("false")
}

/// A value without the blanks Perl drops around a path: its own word split of
/// a command holding no shell metacharacter (`samtools faidx
/// $human_ancestor_location $region`, `LoF.pm` line 539; the `master` `tabix`
/// commands) drops them around the path, so ` /x.fa` opens `/x.fa`; a value
/// equal to its unpadded self carries none of them. The blank set is Perl's:
/// space, tab, line feed, vertical tab, form feed, carriage return. A line
/// feed is a shell metacharacter to Perl as well: a
/// path padded with one sends the command through the shell, which ends it at
/// the line feed before the path, so Ensembl VEP's plugin opens no file there
/// and reads an empty result (no `ANC_ALLELE`, a GERP total of 0), where
/// vep-rs opens the path.
fn perl_unpadded(value: &str) -> &str {
    value.trim_matches(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r'))
}

/// Parse the per-consequence `exon`/`intron` field, `n/total` or, for a variant
/// overlapping two features, `first-last/total`, into `(n, total)`, both 1-based,
/// with `first` as the rank of a range: LoF.pm numifies the field after its
/// `split /\//` (`$exon_idx--`, `$intron_number--`), which reads the leading
/// integer of `6-7`. Returns `None` if absent or malformed.
fn parse_rank_of_total(s: Option<&str>) -> Option<(u32, u32)> {
    let s = s?;
    let (ranks, total) = s.split_once('/')?;
    let first = ranks.split_once('-').map_or(ranks, |(first, _last)| first);
    Some((first.trim().parse().ok()?, total.trim().parse().ok()?))
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
    /// `feature_seq` ref C -> alt T against a GC intron start: the alleles on the
    /// transcript's strand, so a record given on the reverse strand is read on the
    /// forward strand first.
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
        let (forward_ref, forward_alt) = variant.forward_strand_alleles();
        let (mut r, mut a) = (
            forward_ref.first().copied()?.to_ascii_uppercase(),
            forward_alt.first().copied()?.to_ascii_uppercase(),
        );
        if ctx.strand == -1 {
            r = complement(r);
            a = complement(a);
        }
        Some(r == b'C' && a == b'T')
    }

    /// The GERP source Ensembl VEP's LoF plugin reads when none is named: on
    /// GRCh37 master's fixed `gerp_file` path (master `LoF.pm` line 132);
    /// otherwise the bigWig under `loftee_path` by its upstream name (`LoF.pm`
    /// line 130 appends the name to the value as given, so a `loftee_path`
    /// without its trailing slash names a file that does not exist).
    fn default_gerp_source(&self) -> GerpSource {
        match self.assembly {
            Some(vep_core::assembly::Assembly::Grch37) => {
                GerpSource::TabixTsv(GRCH37_DEFAULT_GERP_FILE.to_string())
            }
            _ => {
                let loftee_path = self.loftee_path.as_deref().unwrap_or(DEFAULT_LOFTEE_PATH);
                GerpSource::BigWig(format!("{loftee_path}{GRCH38_DEFAULT_GERP_BIGWIG}"))
            }
        }
    }

    /// The one GERP key Ensembl VEP's LoF plugin reads for the assembly:
    /// `gerp_bigwig` on the `grch38` branch (`LoF.pm` line 246), `gerp_file` on
    /// `master` (master `LoF.pm` line 141 builds the `tabix` command line 258
    /// runs from it). Its file wins over every other key's.
    fn perl_gerp_key(&self) -> &'static str {
        match self.assembly {
            Some(vep_core::assembly::Assembly::Grch37) => "gerp_file",
            _ => "gerp_bigwig",
        }
    }
}

/// The GERP source opened once at `init` and read on every stop-gained or
/// frameshift row, so a broken path fails at init instead of silently
/// disabling END_TRUNC for the whole run, and the tabix index and the bigWig
/// header are read once per run.
struct GerpReader {
    data: GerpData,
    /// The file `data` reads, named when a read of it fails a row.
    path: String,
    /// Set by the first interval on a chromosome the file lacks, which counts
    /// 0 as Ensembl VEP's plugin counts it, so the run warns once, not once
    /// per row.
    no_chromosome_warned: AtomicBool,
}

/// The GERP source's open file, per format.
enum GerpData {
    /// `BigWigRead::get_interval` takes `&mut self` and the plugin is queried
    /// through `&self`. The `vep` pipeline runs one batch's plugin phase at a
    /// time, so the lock is uncontended there; batches annotated on several
    /// threads serialise their GERP reads on it.
    BigWig(Mutex<BigWigRead<ReopenableFile>>),
    /// Queried through `&self`, as every tabix-backed plugin queries its file.
    TabixTsv(TabixAnnotator),
}

impl GerpReader {
    fn open(
        source: &GerpSource,
        assembly: Option<vep_core::assembly::Assembly>,
    ) -> Result<Self, PluginError> {
        let word = source.path();
        let path = source.opened_path();
        if !std::path::Path::new(path).exists() {
            return Err(PluginError::Init(format!(
                "LoFTEE: GERP data file not found: '{word}'"
            )));
        }
        let data = match source {
            GerpSource::BigWig(_) => open_bigwig(path)
                .map(|reader| GerpData::BigWig(Mutex::new(reader)))
                .map_err(|e| {
                    PluginError::Init(format!("LoFTEE: failed to open GERP bigWig '{word}': {e}"))
                })?,
            GerpSource::TabixTsv(_) => {
                // Without its tabix index every interval query returns nothing.
                let config = TabixConfig {
                    file_path: path.into(),
                    chr_col: GERP_TSV_CHR_COL,
                    start_col: GERP_TSV_POS_COL,
                    end_col: None,
                    ref_col: None,
                    alt_col: None,
                    zero_based: false,
                };
                TabixAnnotator::open(config)
                    .map(|table| GerpData::TabixTsv(table.on_assembly(assembly)))
                    .map_err(|e| {
                        PluginError::Init(format!(
                            "LoFTEE: failed to open GERP tabix TSV '{word}': {e}. \
                             Expected the upstream GRCh37 GERP_scores.final.sorted.txt.gz \
                             with a sibling .tbi."
                        ))
                    })?
            }
        };
        Ok(Self {
            data,
            path: path.to_string(),
            no_chromosome_warned: AtomicBool::new(false),
        })
    }

    /// The format's name in a warning: `bigWig` or `tabix TSV`.
    fn format_name(&self) -> &'static str {
        match self.data {
            GerpData::BigWig(_) => "bigWig",
            GerpData::TabixTsv(_) => "tabix TSV",
        }
    }

    /// Warns, once per run, that the GERP file holds no chromosome `chr`, so
    /// the interval `chr:lo-hi` counts 0 toward `GERP_DIST`.
    fn warn_no_chromosome_once(&self, chr: &str, lo: u64, hi: u64) {
        if !self.no_chromosome_warned.swap(true, Ordering::Relaxed) {
            let format = self.format_name();
            warn!(
                "LoFTEE: the GERP {format} {} holds no chromosome {chr}; {chr}:{lo}-{hi} counts \
                 0 toward GERP_DIST, as Ensembl VEP's plugin counts it, and later intervals on a \
                 chromosome the {format} lacks are not reported in this run",
                self.path
            );
        }
    }
}

/// Why one interval of the GERP walk has no total, by what Ensembl VEP's
/// plugin does in its place.
#[derive(Debug)]
enum GerpReadFailure {
    /// The file lacks the chromosome: the interval counts 0 and the plugin
    /// says nothing, as `gerp_dist.pl` gets no feature from a bigWig for one
    /// (`Bio::DB::BigWig::get_seq_stream` returns before any read when
    /// `chromSize` is 0) and no line from `tabix`, which prints nothing for a
    /// region on a sequence its index does not name, so
    /// `get_interval_gerp_tabix` sums 0.
    NoChromosome,
    /// The file could not be read: the row loses its LoF fields under the
    /// plugin's `went wrong` warning, as it does when `Bio::DB::BigWig->new`
    /// dies on a path it cannot open. On master the `tabix` subprocess merges
    /// a read error into the lines `get_interval_gerp_tabix` sums, so the
    /// message's third field reads as 0 there: a defect, not a count, and the
    /// row fails here too. Only an error the bigWig library returns arrives
    /// here; a compressed block it cannot inflate ends the run inside the
    /// library, as the same block ends Ensembl VEP's run.
    FailsRow(String),
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

    /// No LoFTEE word names a remote data file: a word holding a URL holds the
    /// colon of its `://`, so `LoF.pm` splits it into three or more fields and
    /// drops it (lines 68-69), or into a two-field word whose key is the scheme
    /// or holds `=`, which it stores and never reads; `init` reads each word the
    /// same way and logs the ones it drops, so the plugin runs without the word,
    /// as Ensembl VEP's does.
    fn remote_data_file<'a>(&self, _param: &'a str) -> Option<&'a str> {
        None
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        // `LoF.pm` keeps one hash entry per key, its last value, and reads the
        // one key its branch names, so each GERP key holds its own candidate
        // here, in the order the keys were last named.
        let mut gerp_candidates: IndexMap<&str, GerpSource> = IndexMap::new();
        // The two keys that ask for what vep-rs lacks are judged by that last
        // value too, once every word is read.
        let mut conservation_file: Option<&str> = None;
        let mut run_splice_predictions: Option<&str> = None;
        let grch37 = matches!(self.assembly, Some(vep_core::assembly::Assembly::Grch37));
        for param in params {
            let Some((key, value)) = split_parameter(param) else {
                warn_unread_parameters(
                    "LoFTEE",
                    std::slice::from_ref(param),
                    "a key:value word of exactly two colon-separated fields, as Ensembl VEP's \
                     LoF plugin does, or a key=value word",
                );
                continue;
            };
            match key {
                // master's `LoF.pm`, the plugin of a GRCh37 run, reads none of the
                // three: `use_gerp_end_trunc` stands in one default that never
                // runs (line 79) and is tested nowhere, the percentile test that
                // read `filter_position` is commented out (lines 74 and 253), and
                // END_TRUNC compares the GERP total with the literal 180 (line 268).
                "use_gerp_end_trunc" | "filter_position" | "gerp_end_trunc_cutoff" if grch37 => {
                    warn!(
                        "LoFTEE: ignoring '{key}' on GRCh37: Ensembl VEP's LoF plugin for that \
                         assembly never reads it and judges END_TRUNC by the GERP rule at its \
                         own cutoff of 180"
                    );
                }
                // master's `LoF.pm` writes `LoF_info` while this key is false by Perl's
                // string truth (lines 243 and 339); the grch38 branch has no such key.
                "skip_lof_info" if grch37 => {
                    self.skip_lof_info = perl_true(value);
                }
                // master's `LoF.pm` names its `tabix` program with this key (lines 130,
                // 138 and 141); the GERP TSV is read in-process here, so the value has
                // no use. The grch38 branch has no such key.
                "tabix_path" if grch37 => {}
                // `LoF.pm` stores these three values as the text after the colon
                // and reads each in one numeric comparison (lines 297, 241 and
                // 254), so each is Perl's number of its text, blanks around it
                // skipped as Perl skips them, and no shape is refused.
                "min_intron_size" => {
                    self.min_intron_size = perl_numeric("LoFTEE", key, value);
                }
                "filter_position" => {
                    self.filter_position = perl_numeric("LoFTEE", key, value);
                }
                "gerp_end_trunc_cutoff" => {
                    self.gerp_end_trunc_cutoff = perl_numeric("LoFTEE", key, value);
                    self.gerp_cutoff_explicit = true;
                }
                "use_gerp_end_trunc" => {
                    self.use_gerp_end_trunc =
                        !(value == "0" || value.eq_ignore_ascii_case("false"));
                }
                "check_complete_cds" => {
                    // LoF.pm:268 `lc(...) eq 'true'`: any case of `true` is on, `1` is off.
                    self.check_complete_cds = value.eq_ignore_ascii_case("true");
                }
                // `gerp_file` is master's own key for the same tabix TSV.
                "gerp_bigwig" | "gerp_tabix" | "gerp_scores" | "gerp_file" | "gerp" => {
                    gerp_candidates.shift_remove(key);
                    if !is_off_word(value) {
                        let source = match key {
                            "gerp_bigwig" => GerpSource::BigWig(value.to_string()),
                            "gerp" => GerpSource::from_path(value),
                            _ => GerpSource::TabixTsv(value.to_string()),
                        };
                        gerp_candidates.insert(key, source);
                    }
                }
                "human_ancestor_fa" => {
                    // LoF.pm:312 tests the value whole, so a padded ` false` is a path.
                    self.human_ancestor_fa = (!is_off_word(value)).then(|| value.to_string());
                }
                "conservation_file" => conservation_file = Some(value),
                "run_splice_predictions" => run_splice_predictions = Some(value),
                "loftee_path" => {
                    self.loftee_path = (!value.is_empty()).then(|| value.to_string());
                }
                other if param.contains(':') => {
                    warn!(
                        "LoFTEE: ignoring unknown parameter '{other}' in '{param}': a word \
                         holding ':' is read as key:value at that colon, as Ensembl VEP's LoF \
                         plugin reads it"
                    );
                }
                other => {
                    warn!("LoFTEE: ignoring unknown parameter '{other}'");
                }
            }
        }

        // LoF.pm:274 runs the PhyloCSF lookup under any value but `false` in any
        // letter case, tested whole, so a `false` beside a blank asks for it.
        if let Some(value) = conservation_file.filter(|value| !is_off_word(value)) {
            let padded = if is_off_word(perl_unpadded(value)) {
                format!("; the padded value '{value}' names a database in Ensembl VEP's plugin")
            } else {
                String::new()
            };
            return Err(PluginError::Init(format!(
                "LoFTEE: conservation_file={value} asks for the PhyloCSF items of LoF_flags and LoF_info, which vep-rs does not compute{padded}"
            )));
        }
        // Perl LoFTEE reads the value by string truth (LoF.pm:172 on grch38, :186
        // on master), so only `0` is off there. An empty value never reaches that
        // test: `new` drops a parameter whose `split /:/` yields one field, and the
        // branch default applies, 0 on grch38 and 1 on master (LoF.pm:98), so a
        // GRCh37 run holding no value for the key asks for the layer too.
        match run_splice_predictions {
            Some(value) if perl_true(value) => {
                let padded = if perl_unpadded(value) == value {
                    String::new()
                } else {
                    format!("; the padded value '{value}' is true in Ensembl VEP's plugin")
                };
                return Err(PluginError::Init(format!(
                    "LoFTEE: run_splice_predictions={value} asks for the splice-prediction layer, which vep-rs does not implement{padded}"
                )));
            }
            Some("") | None if grch37 => {
                return Err(PluginError::Init(
                    "LoFTEE: this GRCh37 run asks for the splice-prediction layer, which vep-rs \
                     does not implement: Ensembl VEP's LoF plugin runs it on that assembly unless \
                     run_splice_predictions is 0; pass run_splice_predictions:0"
                        .to_string(),
                ));
            }
            _ => {}
        }

        // The key Ensembl VEP's plugin reads wins; a file named under any other
        // key, the keys the plugin never reads, stands in when that key holds
        // none, the last named first.
        self.gerp_source = gerp_candidates
            .shift_remove(self.perl_gerp_key())
            .or_else(|| gerp_candidates.pop().map(|(_, source)| source));

        if !self.gerp_cutoff_explicit {
            if let Some(assembly) = self.assembly {
                self.gerp_end_trunc_cutoff = match assembly {
                    vep_core::assembly::Assembly::Grch37 => MASTER_GERP_END_TRUNC_CUTOFF,
                    vep_core::assembly::Assembly::Grch38 => DEFAULT_GERP_END_TRUNC_CUTOFF,
                };
            }
        }

        if self.use_gerp_end_trunc {
            let reader = match self.gerp_source {
                Some(ref source) => GerpReader::open(source, self.assembly)?,
                None => {
                    let derived = self.default_gerp_source();
                    let note = match derived {
                        GerpSource::BigWig(_) => DERIVED_GERP_BIGWIG_NOTE,
                        GerpSource::TabixTsv(_) => DERIVED_GERP_TABIX_NOTE,
                    };
                    let reader =
                        GerpReader::open(&derived, self.assembly).map_err(|e| e.with_note(note))?;
                    self.gerp_source = Some(derived);
                    reader
                }
            };
            self.gerp_reader = Some(reader);
        }

        if let Some(ref word) = self.human_ancestor_fa {
            // Kept as given so the refusal quotes the word; the file is opened by the
            // path LoF.pm:539's `samtools faidx` command splits out of it (`perl_unpadded`).
            let path = perl_unpadded(word);
            let fasta = vep_fasta::IndexedFasta::from_path(path).map_err(|e| {
                PluginError::Init(format!(
                    "LoFTEE: failed to open human_ancestor_fa '{word}': {e}. \
                     Expected an indexed FASTA: a plain .fa with its .fai, or a bgzipped \
                     .fa.gz with its .fai and .gzi."
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
        if tc.biotype() != Some("protein_coding") {
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
        let mut filters: Vec<&str> = Vec::new();
        let mut flags: Vec<&str> = Vec::new();
        let mut info: Vec<String> = Vec::new();

        let exon_no = parse_rank_of_total(tc.exon.as_deref());
        let intron_no = parse_rank_of_total(tc.intron.as_deref());

        if other_lof && ctx.cdna_coding_start.is_some() && ctx.cdna_coding_end.is_some() {
            if let Some(pct) = cds_percentile(tc, variant, ctx) {
                // LoF.pm joins the numbers into `LoF_info` in Perl's string
                // context, fifteen significant digits (`PERCENTILE:` . $pct).
                info.push(format!("PERCENTILE:{}", perl_number(pct)));
                if !self.use_gerp_end_trunc {
                    if pct >= 1.0 - self.filter_position {
                        filters.push("END_TRUNC");
                    }
                } else if let Some(ref gerp) = self.gerp_reader {
                    if let Some((gerp_dist, bp_dist)) =
                        gerp_weighted_dist(ctx, chr, variant.start, gerp)?
                    {
                        // master's `LoF.pm`, the plugin Ensembl VEP runs on GRCh37,
                        // pushes `PERCENTILE` after `GERP_DIST` and `BP_DIST` (lines
                        // 259-261); the grch38 branch pushes it first (line 239).
                        let percentile = match self.assembly {
                            Some(vep_core::assembly::Assembly::Grch37) => info.pop(),
                            _ => None,
                        };
                        info.push(format!("GERP_DIST:{}", perl_number(gerp_dist)));
                        info.push(format!("BP_DIST:{bp_dist}"));
                        info.extend(percentile);
                        if let Some((rank, _)) = exon_no {
                            let last_exon_len = last_exon_coding_length(ctx, rank);
                            let d = bp_dist as i64 - last_exon_len;
                            info.push(format!("DIST_FROM_LAST_EXON:{d}"));
                            info.push(format!(
                                "50_BP_RULE:{}",
                                if d <= 50 { "FAIL" } else { "PASS" }
                            ));
                            if d <= 50 && gerp_dist <= self.gerp_end_trunc_cutoff {
                                filters.push("END_TRUNC");
                            }
                        } else {
                            flags.push("NO_EXON_NUMBER");
                        }
                    }
                }
            }
        }

        if other_lof {
            // LoF.pm enters its exon checks only with an exon number (`$other_lof &&
            // $tv->exon_number`), so a variant without one gets no exon filter and
            // its call stays HC; the GERP-weighted END_TRUNC rule above is the one
            // place that flags it NO_EXON_NUMBER.
            if let Some((_, total_exons)) = exon_no {
                if total_exons == 1 {
                    flags.push("SINGLE_EXON");
                } else if self.check_complete_cds && (ctx.cds_start_nf || ctx.cds_end_nf) {
                    filters.push("INCOMPLETE_CDS");
                }
            }
        }

        if let Some((intron_rank, _)) = intron_no {
            match ctx.introns.iter().find(|i| i.rank == intron_rank) {
                None => filters.push("EXON_INTRON_UNDEF"),
                Some(intron) => {
                    let intron_size = intron.end.saturating_sub(intron.start) + 1;
                    info.push(format!("INTRON_SIZE:{intron_size}"));
                    if (intron_size as f64) < self.min_intron_size {
                        filters.push("SMALL_INTRON");
                    }
                    if vep_splice_lof {
                        if has_consequence(tc, &[Consequence::SpliceDonorVariant]) {
                            if let Some(true) = self.gc_to_gt_donor(ctx, variant, chr, intron_rank)
                            {
                                filters.push("GC_TO_GT_DONOR");
                            }
                        }
                        if let Some(true) = self.non_canonical_intron(ctx, chr, intron_rank) {
                            flags.push("NON_CAN_SPLICE");
                        }
                        if is_utr_splice(ctx, variant, true) {
                            filters.push("5UTR_SPLICE");
                        }
                        if is_utr_splice(ctx, variant, false) {
                            filters.push("3UTR_SPLICE");
                        }
                    }
                    if has_consequence(tc, &[Consequence::SpliceAcceptorVariant])
                        && self.nagnag_site(variant, chr, ctx)
                    {
                        flags.push("NAGNAG_SITE");
                    }
                }
            }
        }

        if self.ancestor_fasta.is_some() && self.is_ancestral_allele(variant) {
            filters.push("ANC_ALLELE");
        }

        let confidence = if filters.is_empty() { "HC" } else { "LC" };
        out.insert(FIELD_LOF.to_string(), confidence.to_string());
        if !filters.is_empty() {
            out.insert(FIELD_FILTER.to_string(), filters.join(","));
        }
        if !flags.is_empty() {
            out.insert(FIELD_FLAGS.to_string(), flags.join(","));
        }
        if !info.is_empty() && !self.skip_lof_info {
            out.insert(FIELD_INFO.to_string(), info.join(","));
        }
        Ok(out)
    }
}

impl LofteePlugin {
    /// NAGNAG_SITE: acceptor variant whose ±4 genomic context (9 bp,
    /// strand-oriented) carries `AG.AG` at any offset, the unanchored match of
    /// LoF.pm:498 (`$sequence =~ m/AG.AG/`). LoF.pm:497 gates on the length of
    /// the variant's span expanded by four bases each side being nine, so every
    /// one-base reference span is tested (an SNV, a one-base deletion, a
    /// one-base indel) while an insertion (`start = end + 1`, eight bases) and a
    /// longer span are not.
    fn nagnag_site(&self, variant: &InputVariant, chr: &str, ctx: &LofteeContext) -> bool {
        if variant.start != variant.end {
            return false;
        }
        let Some(fasta) = self.reference_fasta.as_ref() else {
            return false;
        };
        let Some(from) = variant.start.checked_sub(4) else {
            return false;
        };
        let Some(seq) = fasta.sequence(chr, from, variant.start + 4) else {
            return false;
        };
        if seq.len() != 9 {
            return false;
        }
        let oriented = if ctx.strand == -1 { revcomp(&seq) } else { seq };
        oriented
            .windows(5)
            .any(|window| window.starts_with(b"AG") && window.ends_with(b"AG"))
    }

    /// ANC_ALLELE: SNP whose alt allele equals the primate ancestral base
    /// (from `human_ancestor.fa`). SNP-only, mirroring LoFTEE's
    /// `check_for_ancestral_allele` (which ignores indels and multi-base alleles).
    ///
    /// LoFTEE compares `variation_feature_seq`, the alt allele as the input record
    /// writes it, with the forward-strand ancestral base (LoF.pm
    /// `check_for_ancestral_allele`), so a record given on the reverse strand is
    /// compared as written, not through its complement.
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
/// and the cDNA coding span from the context. `cds_position` prints `n-m` with
/// the smaller coordinate first, and Ensembl writes an insertion as `start =
/// end + 1`, so its `cds_end` is the smaller number where a deletion's is the
/// larger.
fn cds_percentile(
    tc: &TranscriptConsequence,
    variant: &InputVariant,
    ctx: &LofteeContext,
) -> Option<f64> {
    let cds_len = ctx.cdna_coding_end? as i64 - ctx.cdna_coding_start? as i64 + 1;
    if cds_len <= 0 {
        return None;
    }
    let cds_pos = tc.cds_position.as_deref()?;
    let coordinates = cds_pos
        .split(['-', '/'])
        .filter_map(|p| p.trim().parse::<i64>().ok());
    let end = if variant.variant_class == VariantClass::Insertion {
        coordinates.min()?
    } else {
        coordinates.max()?
    };
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
/// `None` without a coding region. An interval on a chromosome the file lacks
/// counts 0, warned about once per run; an interval the file could not give
/// fails the row ([`GerpReadFailure`]).
fn gerp_weighted_dist(
    ctx: &LofteeContext,
    chr: &str,
    pos: u64,
    reader: &GerpReader,
) -> Result<Option<(f64, u64)>, PluginError> {
    let (Some(coding_start), Some(coding_end)) = (ctx.coding_region_start, ctx.coding_region_end)
    else {
        return Ok(None);
    };
    let stop_codon_pos = if ctx.strand == 1 {
        coding_end
    } else {
        coding_start
    };

    let mut summer = GerpSummer::of(reader);

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
        match summer.sum_interval(chr, lo, hi) {
            Ok(total) => weighted_dist += total,
            Err(GerpReadFailure::NoChromosome) => reader.warn_no_chromosome_once(chr, lo, hi),
            // The text names the file and the error, never the interval, so a
            // file whose reads fail with one text warns once (the sink writes
            // each text once); a tabix error naming a block offset is its own
            // text.
            Err(GerpReadFailure::FailsRow(reason)) => {
                return Err(PluginError::Run(format!("LoFTEE: {reason}")));
            }
        }
        bp_dist += hi.abs_diff(lo);
    }

    Ok(Some((weighted_dist, bp_dist)))
}

/// Per-format GERP interval summer.
///
/// Both upstream branches walk the same exon intervals; they differ only in how
/// they total one interval:
///
/// - **bigWig** (grch38): `gerp_dist.pl` reads one summary feature over the
///   range and adds `$c->length * binMean($c->score)`, the range's length times
///   the mean of its covered bases, so the total is evaluated here in that order
///   too (`length * (sum / covered)`), which is the same double as the plain
///   sum of the values wherever the range is fully covered and the same text in
///   `LoF_info` even where the two roundings differ. Coordinates are 0-based,
///   half-open `[start, end)`.
/// - **tabix TSV** (master/GRCh37): the file holds one row per base, and
///   `gerp_dist.pl` does `$gerp += $res[2]`, a plain sum of per-base values with
///   no length multiplier. Coordinates are 1-based inclusive.
///
/// One summer serves one `gerp_weighted_dist` call, every interval of its exon
/// walk, over the plugin's [`GerpReader`]: the bigWig reader locked for the
/// call, the tabix annotator borrowed.
enum GerpSummer<'a> {
    BigWig {
        reader: MutexGuard<'a, BigWigRead<ReopenableFile>>,
        path: &'a str,
    },
    TabixTsv(&'a TabixAnnotator),
}

impl<'a> GerpSummer<'a> {
    fn of(reader: &'a GerpReader) -> Self {
        match &reader.data {
            GerpData::BigWig(bigwig) => GerpSummer::BigWig {
                reader: bigwig.lock().unwrap_or_else(PoisonError::into_inner),
                path: &reader.path,
            },
            GerpData::TabixTsv(annotator) => GerpSummer::TabixTsv(annotator),
        }
    }

    /// Total GERP over the 1-based inclusive range `[lo, hi]`.
    fn sum_interval(&mut self, chr: &str, lo: u64, hi: u64) -> Result<f64, GerpReadFailure> {
        match self {
            GerpSummer::BigWig { reader, path } => {
                // bigWig is 0-based half-open, so shift the start and use `hi`
                // directly as the exclusive end.
                let q_start = lo.saturating_sub(1) as u32;
                let q_end = hi as u32;
                let mut sum = 0.0_f64;
                let mut covered = 0.0_f64;
                let values = reader
                    .get_interval(chr, q_start, q_end)
                    .map_err(|e| classify_bigwig_failure(path, e))?;
                for value in values {
                    let value = value.map_err(|e| classify_bigwig_failure(path, e))?;
                    let clipped_start = value.start.max(q_start);
                    let clipped_end = value.end.min(q_end);
                    if clipped_end > clipped_start {
                        let bases = (clipped_end - clipped_start) as f64;
                        sum += bases * value.value as f64;
                        covered += bases;
                    }
                }
                if covered == 0.0 {
                    return Ok(0.0);
                }
                Ok((q_end - q_start) as f64 * (sum / covered))
            }
            GerpSummer::TabixTsv(annotator) => {
                // A chromosome or span the index lacks is `Ok` and empty; an
                // `Err` is a read failure, and its text names the file. The
                // index is asked about the chromosome on an empty span alone,
                // so a span with records pays no lookup.
                let records = annotator.query(chr, lo, hi).map_err(|e| {
                    GerpReadFailure::FailsRow(match e {
                        PluginError::Run(reason) => reason,
                        other => other.to_string(),
                    })
                })?;
                if records.is_empty() && !annotator.holds_chromosome(chr) {
                    return Err(GerpReadFailure::NoChromosome);
                }
                let mut total = 0.0_f64;
                for record in records {
                    if let Some(field) = record.columns.get(GERP_TSV_SCORE_COL) {
                        if let Ok(v) = field.trim().parse::<f64>() {
                            total += v;
                        }
                    }
                }
                Ok(total)
            }
        }
    }
}

/// A chromosome the bigWig lacks counts 0, as `Bio::DB::BigWig` gives
/// `gerp_dist.pl` no feature for one; any error the library returns fails the
/// row. A compressed block the library cannot inflate never reaches here: it
/// ends the run inside the library, as the same block ends Ensembl VEP's run.
fn classify_bigwig_failure(path: &str, error: BBIReadError) -> GerpReadFailure {
    match error {
        BBIReadError::InvalidChromosome(_) => GerpReadFailure::NoChromosome,
        _ => GerpReadFailure::FailsRow(format!("reading {path}: {error}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::ConsequenceList;
    use vep_core::coordinate::Strand;
    use vep_core::transcript::Transcript;

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

    /// A variant overlapping two exons or introns carries `first-last/total`
    /// (BaseTranscriptVariation.pm:700 and :748), which LoF.pm:565-566 and
    /// :511-512 numify to the first feature's rank (`6-7` reads 6) and
    /// LoF.pm:457-458 reads the total of, so the range is a rank, not a
    /// missing number.
    #[test]
    fn parse_rank_of_total_reads_the_first_feature_of_a_range() {
        assert_eq!(parse_rank_of_total(Some("6-7/9")), Some((6, 9)));
        assert_eq!(parse_rank_of_total(Some("6-7/12")), Some((6, 12)));
        assert_eq!(parse_rank_of_total(Some("1-2/2")), Some((1, 2)));
    }

    /// LoF.pm:268 tests `lc($self->{check_complete_cds}) eq 'true'`, so every
    /// letter case of `true` turns the `INCOMPLETE_CDS` filter on and any other
    /// word, `1` included, leaves it off.
    #[test]
    fn check_complete_cds_reads_any_case_of_true_only() {
        for (value, on) in [
            ("true", true),
            ("True", true),
            ("TRUE", true),
            ("tRuE", true),
            ("1", false),
            ("yes", false),
            ("false", false),
        ] {
            let mut p = LofteePlugin::new();
            p.init(&[
                format!("check_complete_cds:{value}"),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
            assert_eq!(p.check_complete_cds, on, "check_complete_cds:{value}");
        }
    }

    /// LoF.pm:312 tests `lc($self->{human_ancestor_fa}) ne 'false'`, so every
    /// letter case of `false` switches the ancestral-allele filter off, and
    /// `init` opens no file. A word that merely starts with `false` is a path
    /// there, read by `samtools faidx`, so it stays a path here and the missing
    /// file is refused at `init`.
    #[test]
    fn human_ancestor_fa_reads_any_case_of_false_as_off() {
        for value in ["false", "False", "FALSE", "fAlSe"] {
            let mut p = LofteePlugin::new();
            p.init(&[
                format!("human_ancestor_fa:{value}"),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap_or_else(|e| panic!("human_ancestor_fa:{value}: {e}"));
            assert_eq!(p.human_ancestor_fa, None, "human_ancestor_fa:{value}");
            assert!(p.ancestor_fasta.is_none(), "human_ancestor_fa:{value}");
        }

        let mut path = LofteePlugin::new();
        let err = path
            .init(&[
                "human_ancestor_fa:false.fa".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .expect_err("a path starting with `false` is a path")
            .to_string();
        assert!(
            err.contains("failed to open human_ancestor_fa 'false.fa'"),
            "unexpected error: {err}"
        );
    }

    /// LoF.pm:274 tests `lc($self->{conservation_file}) ne 'false'` before the
    /// PhyloCSF lookup, so every letter case of `false` leaves the lookup off,
    /// as the plugin's default does (line 119), and so does the empty value of
    /// a `key=` word: a stop-gained row with an exon number, the row the lookup
    /// runs on (line 262), is the one a run without the parameter writes.
    #[test]
    fn conservation_file_reads_any_case_of_false_as_off() {
        let row_under = |word: Option<&str>| {
            let mut params = vec!["use_gerp_end_trunc=0".to_string()];
            params.extend(word.map(str::to_string));
            let mut p = LofteePlugin::new();
            p.init(&params).unwrap_or_else(|e| panic!("{word:?}: {e}"));
            let (tc, v) = gerp_row("21");
            p.run(&tc, &v).unwrap()
        };
        let without = row_under(None);
        assert_eq!(
            without
                .iter()
                .map(|(field, value)| (field.as_str(), value.as_str()))
                .collect::<Vec<_>>(),
            [
                ("LoF", "HC"),
                ("LoF_flags", "SINGLE_EXON"),
                ("LoF_info", "PERCENTILE:0.166666666666667"),
            ]
        );
        for word in [
            "conservation_file:false",
            "conservation_file:False",
            "conservation_file:FALSE",
            "conservation_file=false",
            "conservation_file=",
        ] {
            assert_eq!(row_under(Some(word)), without, "{word}");
        }
    }

    /// Under any other `conservation_file` LoF.pm writes the PhyloCSF items: the
    /// `LoF_info` items `ANN_ORF` and `MAX_ORF`, or `PHYLOCSF_TOO_SHORT` for an
    /// exon its database lacks, and the `LoF_flags` `PHYLOCSF_WEAK` or
    /// `PHYLOCSF_UNLIKELY_ORF` (lines 274-285). vep-rs computes none of them, so
    /// a value naming a database, a file or the word `mysql` (line 122), is
    /// refused at initialization in either grammar, quoted as given.
    #[test]
    fn conservation_file_naming_a_database_is_refused() {
        for (param, value) in [
            (
                "conservation_file:/data/phylocsf.sqlite",
                "/data/phylocsf.sqlite",
            ),
            (
                "conservation_file=/data/phylocsf.sqlite",
                "/data/phylocsf.sqlite",
            ),
            (
                "conservation_file: /data/phylocsf.sqlite",
                " /data/phylocsf.sqlite",
            ),
            ("conservation_file:mysql", "mysql"),
        ] {
            let err = LofteePlugin::new()
                .init(&[param.to_string(), "use_gerp_end_trunc=0".to_string()])
                .expect_err(param);
            match err {
                PluginError::Init(message) => assert_eq!(
                    message,
                    format!(
                        "LoFTEE: conservation_file={value} asks for the PhyloCSF items of \
                         LoF_flags and LoF_info, which vep-rs does not compute"
                    ),
                    "{param}"
                ),
                other => panic!("{param}: expected PluginError::Init, got {other:?}"),
            }
        }
    }

    /// `LoF.pm` keeps a value as `split /:/` yields it (lines 68-70) and tests
    /// it whole (lines 121 and 274), so a `false` beside a blank is no off word
    /// there: it names the database line 126 opens and line 275 reads, and the
    /// word is refused as one.
    #[test]
    fn conservation_file_with_a_blank_beside_its_false_names_a_database_as_in_perl() {
        for param in ["conservation_file: false", "conservation_file:False\t"] {
            let err = LofteePlugin::new()
                .init(&[param.to_string(), "use_gerp_end_trunc=0".to_string()])
                .expect_err(param);
            let value = param.split_once(':').unwrap().1;
            match err {
                PluginError::Init(message) => assert_eq!(
                    message,
                    format!(
                        "LoFTEE: conservation_file={value} asks for the PhyloCSF items of \
                         LoF_flags and LoF_info, which vep-rs does not compute; the padded \
                         value '{value}' names a database in Ensembl VEP's plugin"
                    ),
                    "{param}"
                ),
                other => panic!("{param}: expected PluginError::Init, got {other:?}"),
            }
        }
    }

    /// vep-rs's off words for `use_gerp_end_trunc` are `0` and `false` in any
    /// letter case; every other value, an empty one included, keeps the GERP
    /// rule, the plugin's default.
    #[test]
    fn use_gerp_end_trunc_reads_any_case_of_false_as_off() {
        for (value, on) in [
            ("0", false),
            ("false", false),
            ("False", false),
            ("FALSE", false),
            ("fAlSe", false),
            ("1", true),
            ("true", true),
            ("no", true),
            ("", true),
        ] {
            let mut p = LofteePlugin::new();
            p.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
            let result = p.init(&[
                format!("use_gerp_end_trunc={value}"),
                "gerp_bigwig=/nonexistent/gerp.bw".to_string(),
            ]);
            assert_eq!(p.use_gerp_end_trunc, on, "use_gerp_end_trunc={value}");
            assert_eq!(result.is_err(), on, "use_gerp_end_trunc={value}");
        }
    }

    /// The GERP-source keys read the same off word as the other data-file
    /// parameters: `false` in any letter case, or an empty value, names no
    /// file, under both grammars, so with the GERP rule off `init` opens nothing
    /// and keeps no source, and with it on the default source is derived; a
    /// word that is not the off word stays a path, a padded ` false` included,
    /// as it is the file name `LoF.pm` stores: the source holds the word as
    /// given and opens the bigWig by it, the tabix path without the blanks
    /// `master`'s `tabix` command split drops. `LoF.pm` keeps one hash entry
    /// per key, so an off word on one key leaves a source another key named.
    #[test]
    fn gerp_source_keys_read_any_case_of_false_as_no_file() {
        const KEYS: [&str; 5] = [
            "gerp_bigwig",
            "gerp_tabix",
            "gerp_scores",
            "gerp_file",
            "gerp",
        ];
        for key in KEYS {
            for word in ["false", "False", "FALSE", "fAlSe", ""] {
                for sep in ['=', ':'] {
                    let mut p = LofteePlugin::new();
                    p.init(&[
                        format!("{key}{sep}{word}"),
                        "use_gerp_end_trunc=0".to_string(),
                    ])
                    .unwrap_or_else(|e| panic!("{key}{sep}{word}: {e}"));
                    assert!(p.gerp_source.is_none(), "{key}{sep}{word}");
                }
            }
            for path_word in ["/data/false.bw", " false", " /data/false.bw\t", " "] {
                let mut path = LofteePlugin::new();
                path.init(&[
                    format!("{key}={path_word}"),
                    "use_gerp_end_trunc=0".to_string(),
                ])
                .unwrap();
                let bigwig = key == "gerp_bigwig"
                    || (key == "gerp" && perl_unpadded(path_word).ends_with(".bw"));
                assert_eq!(
                    matches!(path.gerp_source, Some(GerpSource::BigWig(_))),
                    bigwig,
                    "{key}={path_word}: {:?}",
                    path.gerp_source
                );
                assert_eq!(
                    path.gerp_source.as_ref().map(GerpSource::path),
                    Some(path_word),
                    "{key}={path_word}"
                );
                let opened = if bigwig {
                    path_word
                } else {
                    perl_unpadded(path_word)
                };
                assert_eq!(
                    path.gerp_source.as_ref().map(GerpSource::opened_path),
                    Some(opened),
                    "{key}={path_word}"
                );
            }
        }

        let mut same_key = LofteePlugin::new();
        same_key
            .init(&[
                "gerp_bigwig=/data/a.bw".to_string(),
                "gerp_bigwig:FALSE".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        assert!(same_key.gerp_source.is_none());
        let mut other_key = LofteePlugin::new();
        other_key
            .init(&[
                "gerp_bigwig=/data/a.bw".to_string(),
                "gerp_file:FALSE".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        assert_eq!(
            other_key.gerp_source.as_ref().map(GerpSource::path),
            Some("/data/a.bw")
        );

        let dir = tempfile::tempdir().unwrap();
        let bigwig = dir.path().join(GRCH38_DEFAULT_GERP_BIGWIG);
        std::fs::copy(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
                "../../tests/golden/116/GRCh38-plugin-loftee/plugin_data/\
                 gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            ),
            &bigwig,
        )
        .unwrap();
        let mut derived = LofteePlugin::new();
        derived
            .init(&[
                format!("loftee_path:{}/", dir.path().display()),
                "gerp_bigwig:False".to_string(),
            ])
            .unwrap();
        assert_eq!(
            derived.gerp_source.as_ref().map(GerpSource::path),
            Some(bigwig.to_str().unwrap())
        );
        assert!(derived.gerp_reader.is_some());
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
        assert_eq!(p.min_intron_size, 20.0);
        assert!(!p.use_gerp_end_trunc);
        assert_eq!(p.gerp_end_trunc_cutoff, -60.0);
        assert_eq!(
            p.gerp_source.as_ref().map(|s| s.path()),
            Some("/path/to/loftee/grch38/gerp.bw")
        );
    }

    /// The GRCh37 bundle ships GERP as a tabix per-base TSV, not a bigWig, so
    /// `gerp=` must classify by extension rather than assuming bigWig. The
    /// extension is read inside the blanks around the path, which the source
    /// keeps as given; the path opened loses every blank of Perl's set around
    /// a tabix path, a carriage return included, and none around a bigWig.
    #[test]
    fn gerp_source_classifies_by_extension() {
        let cr_padded_tabix = GerpSource::from_path("/d/x.txt.gz\r");
        assert!(matches!(&cr_padded_tabix, GerpSource::TabixTsv(p) if p == "/d/x.txt.gz\r"));
        assert_eq!(cr_padded_tabix.opened_path(), "/d/x.txt.gz");
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
        let padded_bigwig = GerpSource::from_path(" /d/GERP.BW\t");
        assert!(matches!(&padded_bigwig, GerpSource::BigWig(p) if p == " /d/GERP.BW\t"));
        assert_eq!(padded_bigwig.opened_path(), " /d/GERP.BW\t");
        let space_padded_tabix = GerpSource::from_path(" /d/x.txt.gz ");
        assert!(matches!(&space_padded_tabix, GerpSource::TabixTsv(p) if p == " /d/x.txt.gz "));
        assert_eq!(space_padded_tabix.opened_path(), "/d/x.txt.gz");
    }

    /// The word a GRCh37 run passes: master's `LoF.pm` runs the
    /// splice-prediction layer unless the key holds `0` (lines 98 and 186).
    const SPLICE_PREDICTIONS_OFF: &str = "run_splice_predictions:0";

    /// A one-row tabix GERP TSV of the layout the GRCh37 bundle ships, under `dir`.
    fn fixture_gerp_tabix(dir: &std::path::Path) -> std::path::PathBuf {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        write_tabix_fixture(
            dir,
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", "21\t1500\t1.5"],
            IndexSpec::standard(),
        )
    }

    /// `gerp_tabix=` selects the per-base TSV reader, the format the GRCh37
    /// bundle ships, and a GRCh37 run opens it at `init`.
    #[test]
    fn accepts_a_tabix_gerp_source_for_grch37() {
        let dir = tempfile::tempdir().unwrap();
        let mut p = LofteePlugin::new();
        p.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        p.init(&[
            format!("gerp_tabix={}", fixture_gerp_tabix(dir.path()).display()),
            SPLICE_PREDICTIONS_OFF.to_string(),
        ])
        .unwrap();
        assert!(matches!(p.gerp_source, Some(GerpSource::TabixTsv(_))));
        assert!(p.gerp_reader.is_some());
    }

    /// The END_TRUNC cutoff differs in sign between the branches (+180 master vs
    /// -58 grch38). Applying the grch38 value to GRCh37 data would make END_TRUNC
    /// effectively never fire, so it is derived from the assembly.
    #[test]
    fn end_trunc_cutoff_follows_the_assembly() {
        let percentile_rule = ["use_gerp_end_trunc=0".to_string()];
        let dir = tempfile::tempdir().unwrap();
        let mut g37 = LofteePlugin::new();
        g37.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        g37.init(&[
            format!("gerp_tabix={}", fixture_gerp_tabix(dir.path()).display()),
            SPLICE_PREDICTIONS_OFF.to_string(),
        ])
        .unwrap();
        assert_eq!(g37.gerp_end_trunc_cutoff, MASTER_GERP_END_TRUNC_CUTOFF);

        let mut g38 = LofteePlugin::new();
        g38.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        g38.init(&percentile_rule).unwrap();
        assert_eq!(g38.gerp_end_trunc_cutoff, DEFAULT_GERP_END_TRUNC_CUTOFF);

        let mut none = LofteePlugin::new();
        none.init(&percentile_rule).unwrap();
        assert_eq!(none.gerp_end_trunc_cutoff, DEFAULT_GERP_END_TRUNC_CUTOFF);
    }

    /// An explicit cutoff survives assembly derivation on GRCh38, where the
    /// branch reads the key (`LoF.pm` line 254); on GRCh37 master's literal 180
    /// stands (its line 268).
    #[test]
    fn explicit_cutoff_overrides_the_assembly_default() {
        let explicit = [
            "gerp_end_trunc_cutoff=42.5".to_string(),
            "use_gerp_end_trunc=0".to_string(),
        ];
        let mut g38 = LofteePlugin::new();
        g38.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        g38.init(&explicit).unwrap();
        assert_eq!(g38.gerp_end_trunc_cutoff, 42.5);

        let dir = tempfile::tempdir().unwrap();
        let mut g37 = LofteePlugin::new();
        g37.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        g37.init(&[
            "gerp_end_trunc_cutoff=42.5".to_string(),
            format!("gerp_tabix={}", fixture_gerp_tabix(dir.path()).display()),
            SPLICE_PREDICTIONS_OFF.to_string(),
        ])
        .unwrap();
        assert_eq!(g37.gerp_end_trunc_cutoff, MASTER_GERP_END_TRUNC_CUTOFF);
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

    /// With the GERP rule on and no GERP source named, `init` reads the source
    /// Ensembl VEP's LoF plugin reads: `loftee_path` (`/vep/loftee/` when absent,
    /// `LoF.pm` line 83) followed by `gerp_conservation_scores.homo_sapiens.GRCh38.bw`
    /// with no separator (`LoF.pm` line 130), so a `loftee_path` without its
    /// trailing slash names a file that does not exist and is refused with the
    /// note saying where the path came from; on GRCh37 master's fixed
    /// `/vep/loftee/GERP_scores.final.sorted.txt.gz` (master `LoF.pm` line 132).
    /// `use_gerp_end_trunc=0` derives nothing where the key is read; on GRCh37
    /// it is not, so the default is derived under it too. `gerp_file`, master's
    /// own key, names a tabix TSV validated as `gerp_tabix` is.
    #[test]
    fn init_derives_the_default_gerp_source_from_loftee_path() {
        let dir = tempfile::tempdir().unwrap();
        let bigwig = dir.path().join(GRCH38_DEFAULT_GERP_BIGWIG);
        std::fs::copy(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
                "../../tests/golden/116/GRCh38-plugin-loftee/plugin_data/\
                 gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            ),
            &bigwig,
        )
        .unwrap();
        let init_message = |result: Result<(), PluginError>| -> String {
            match result {
                Err(PluginError::Init(message)) => message,
                other => panic!("expected PluginError::Init, got {other:?}"),
            }
        };

        let mut with_slash = LofteePlugin::new();
        with_slash
            .init(&[format!("loftee_path:{}/", dir.path().display())])
            .unwrap();
        assert!(
            matches!(&with_slash.gerp_source, Some(GerpSource::BigWig(p)) if *p == bigwig.display().to_string()),
            "{:?}",
            with_slash.gerp_source
        );

        let without_slash = format!("{}{GRCH38_DEFAULT_GERP_BIGWIG}", dir.path().display());
        assert_eq!(
            init_message(
                LofteePlugin::new().init(&[format!("loftee_path={}", dir.path().display())])
            ),
            format!(
                "LoFTEE: GERP data file not found: '{without_slash}'; {DERIVED_GERP_BIGWIG_NOTE}"
            )
        );

        // The bare word's paths come from the constants alone; their refusals
        // are asserted only where the plugin's own default directory is absent.
        let bare = LofteePlugin::new();
        assert!(matches!(
            bare.default_gerp_source(),
            GerpSource::BigWig(p) if p == "/vep/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
        ));
        let mut g37 = LofteePlugin::new();
        g37.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        assert!(matches!(
            g37.default_gerp_source(),
            GerpSource::TabixTsv(p) if p == "/vep/loftee/GERP_scores.final.sorted.txt.gz"
        ));
        if !std::path::Path::new(DEFAULT_LOFTEE_PATH).exists() {
            assert_eq!(
                init_message(LofteePlugin::new().init(&[])),
                format!(
                    "LoFTEE: GERP data file not found: \
                     '/vep/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw'; \
                     {DERIVED_GERP_BIGWIG_NOTE}"
                )
            );
            for words in [
                &[SPLICE_PREDICTIONS_OFF.to_string()][..],
                &[
                    SPLICE_PREDICTIONS_OFF.to_string(),
                    "use_gerp_end_trunc=0".to_string(),
                ][..],
            ] {
                let mut g37 = LofteePlugin::new();
                g37.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
                assert_eq!(
                    init_message(g37.init(words)),
                    format!(
                        "LoFTEE: GERP data file not found: \
                         '/vep/loftee/GERP_scores.final.sorted.txt.gz'; {DERIVED_GERP_TABIX_NOTE}"
                    ),
                    "{words:?}"
                );
            }
        }

        let mut percentile = LofteePlugin::new();
        percentile
            .init(&["use_gerp_end_trunc=0".to_string()])
            .unwrap();
        assert!(percentile.gerp_source.is_none());

        let mut gerp_file = LofteePlugin::new();
        gerp_file
            .init(&[
                "gerp_file=/d/GERP_scores.final.sorted.txt.gz".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        assert!(matches!(
            gerp_file.gerp_source,
            Some(GerpSource::TabixTsv(ref p)) if p == "/d/GERP_scores.final.sorted.txt.gz"
        ));
        assert_eq!(
            init_message(
                LofteePlugin::new()
                    .init(&["gerp_file=/d/GERP_scores.final.sorted.txt.gz".to_string()])
            ),
            "LoFTEE: GERP data file not found: '/d/GERP_scores.final.sorted.txt.gz'"
        );
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

        let reader = GerpReader::open(
            &GerpSource::TabixTsv(path.to_string_lossy().to_string()),
            None,
        )
        .expect("open tabix GERP");
        let mut summer = GerpSummer::of(&reader);

        // 1.5 + 2.0 + (-0.5) = 3.0, a plain sum with no length weighting.
        let total = summer.sum_interval("21", 100, 102).unwrap();
        assert!(
            (total - 3.0).abs() < 1e-9,
            "expected per-base sum 3.0, got {total}"
        );

        let partial = summer.sum_interval("21", 100, 101).unwrap();
        assert!(
            (partial - 3.5).abs() < 1e-9,
            "expected 3.5 for the first two bases, got {partial}"
        );
    }

    /// `LoF.pm` keeps one hash entry per key, its last value (lines 67-72),
    /// and `run` reads the one key its branch names: `gerp_bigwig` on `grch38`
    /// (line 246), `gerp_file` on `master` (lines 132 and 258); a word on any
    /// other key changes nothing there. Here each key holds its own candidate,
    /// an off word clears that key alone, the assembly's key wins when it holds
    /// a file, and the last-named of the other keys stands in when it does not:
    /// on GRCh38 every order of `gerp_bigwig:A`, `gerp_file:T` and
    /// `gerp_file:false` names A.
    #[test]
    fn gerp_source_keys_resolve_per_key_as_the_assembly_reads_them() {
        use vep_core::assembly::Assembly;
        const A: &str = "gerp_bigwig:/data/a.bw";
        const T: &str = "gerp_file:/data/t.gz";
        const T_OFF: &str = "gerp_file:false";
        let source_of = |assembly: Option<Assembly>, words: &[&str]| -> Option<String> {
            let mut p = LofteePlugin::new();
            p.set_assembly(assembly);
            let mut params: Vec<String> = words.iter().map(|w| w.to_string()).collect();
            params.push("use_gerp_end_trunc=0".to_string());
            if assembly == Some(Assembly::Grch37) {
                params.push(SPLICE_PREDICTIONS_OFF.to_string());
            }
            let result = p.init(&params);
            let source = p.gerp_source.as_ref().map(|s| s.path().to_string());
            match (result, assembly) {
                (Ok(()), _) => {}
                // A GRCh37 run opens its source whatever `use_gerp_end_trunc`
                // says, and these paths name no file: the refusal names the
                // source the keys resolved to.
                (Err(e), Some(Assembly::Grch37)) => {
                    let path = source.as_deref().unwrap_or_default();
                    assert!(
                        e.to_string()
                            .contains(&format!("GERP data file not found: '{path}'")),
                        "{words:?}: {e}"
                    );
                }
                (Err(e), _) => panic!("{words:?}: {e}"),
            }
            source
        };
        let permutations: [[&str; 3]; 6] = [
            [A, T, T_OFF],
            [A, T_OFF, T],
            [T, A, T_OFF],
            [T, T_OFF, A],
            [T_OFF, A, T],
            [T_OFF, T, A],
        ];
        for words in &permutations {
            for assembly in [None, Some(Assembly::Grch38)] {
                assert_eq!(
                    source_of(assembly, words),
                    Some("/data/a.bw".to_string()),
                    "{words:?} {assembly:?}"
                );
            }
            // master reads `gerp_file`, whose own last word decides.
            let file_named_last =
                words.iter().rposition(|w| *w == T) > words.iter().rposition(|w| *w == T_OFF);
            let expected = if file_named_last {
                "/data/t.gz"
            } else {
                "/data/a.bw"
            };
            assert_eq!(
                source_of(Some(Assembly::Grch37), words),
                Some(expected.to_string()),
                "{words:?} GRCh37"
            );
        }

        for words in [[A, T], [T, A]] {
            assert_eq!(
                source_of(Some(Assembly::Grch38), &words),
                Some("/data/a.bw".to_string()),
                "{words:?}"
            );
            assert_eq!(
                source_of(Some(Assembly::Grch37), &words),
                Some("/data/t.gz".to_string()),
                "{words:?}"
            );
        }
        // An off word on the assembly's key leaves the other key's file.
        for assembly in [Some(Assembly::Grch38), Some(Assembly::Grch37)] {
            assert_eq!(
                source_of(assembly, &[T, A, "gerp_bigwig:false"]),
                Some("/data/t.gz".to_string()),
                "{assembly:?}"
            );
        }
        assert_eq!(
            source_of(None, &[A, "gerp_bigwig:false", A]),
            Some("/data/a.bw".to_string())
        );
        // Among the keys Ensembl VEP's plugin never reads, the last named wins.
        const T1: &str = "gerp_tabix:/data/t1.gz";
        const T2: &str = "gerp:/data/t2.bw";
        assert_eq!(source_of(None, &[T1, T2]), Some("/data/t2.bw".to_string()));
        assert_eq!(source_of(None, &[T2, T1]), Some("/data/t1.gz".to_string()));
        assert_eq!(
            source_of(None, &[T1, T2, "gerp:false"]),
            Some("/data/t1.gz".to_string())
        );
        assert_eq!(
            source_of(None, &[T1, T2, "gerp_tabix:FALSE"]),
            Some("/data/t2.bw".to_string())
        );
        assert_eq!(source_of(None, &["gerp_bigwig:false", T_OFF]), None);
    }

    /// The committed GERP bigWig: chromosome `21` alone, the corpus's exons.
    fn fixture_gerp_bigwig() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
            "../../tests/golden/116/GRCh38-plugin-loftee/plugin_data/\
             gerp_conservation_scores.homo_sapiens.GRCh38.bw",
        )
    }

    /// A stop-gained row on the one coding exon `[1000, 2000]` of `coding_ctx`,
    /// the variant at 1500 on `chr`: the GERP walk totals one interval,
    /// `chr:1500-2000`, so `BP_DIST` is 500 and the distance from the last exon
    /// -500.
    fn gerp_row(chr: &str) -> (TranscriptConsequence, InputVariant) {
        let mut ctx = coding_ctx(1);
        ctx.exons = vec![vep_core::transcript::Exon {
            start: 1000,
            end: 2000,
            rank: 1,
            phase: 0,
            end_phase: 0,
        }];
        let tc = TranscriptConsequence {
            transcript: Some(transcript_of_biotype("protein_coding")),
            consequences: smallvec![Consequence::StopGained],
            exon: Some("1/1".to_string()),
            cds_position: Some("150".to_string()),
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        };
        let variant = InputVariant::new(chr.into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        (tc, variant)
    }

    /// A plugin over a copy of the committed GERP bigWig, cut short after `init`
    /// opened it so the index every interval read starts from is gone and each
    /// read fails; with the copy's path and its directory, kept alive with the
    /// plugin.
    fn plugin_over_a_cut_gerp_bigwig() -> (LofteePlugin, std::path::PathBuf, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let bigwig = dir.path().join("gerp.bw");
        std::fs::copy(fixture_gerp_bigwig(), &bigwig).unwrap();
        let mut p = LofteePlugin::new();
        p.init(&[format!("gerp_bigwig={}", bigwig.display())])
            .unwrap();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&bigwig)
            .unwrap()
            .set_len(64)
            .unwrap();
        (p, bigwig, dir)
    }

    /// `gerp_dist.pl` totals an interval through `Bio::DB::BigWig::features`,
    /// which answers a chromosome the bigWig lacks with no feature at all
    /// (`get_seq_stream` returns before any read when `chromSize` is 0), so the
    /// interval counts 0 toward `GERP_DIST` without a word. vep-rs counts the
    /// same 0 and says so once per run, not once per row.
    #[test]
    fn a_chromosome_the_gerp_bigwig_lacks_counts_zero_and_is_warned_once() {
        use crate::test_fixtures::Warnings;
        let mut p = LofteePlugin::new();
        p.init(&[format!("gerp_bigwig={}", fixture_gerp_bigwig().display())])
            .unwrap();
        let warnings = Arc::new(Warnings::default());
        let outputs = tracing::subscriber::with_default(Arc::clone(&warnings), || {
            let (tc, v) = gerp_row("22");
            let first = p.run(&tc, &v).unwrap();
            let (tc, v) = gerp_row("22");
            let second = p.run(&tc, &v).unwrap();
            [first, second]
        });
        for out in &outputs {
            assert_eq!(out.get("LoF").map(String::as_str), Some("HC"), "{out:?}");
            assert_eq!(
                out.get("LoF_info").map(String::as_str),
                Some(
                    "PERCENTILE:0.166666666666667,GERP_DIST:0,BP_DIST:500,\
                     DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                )
            );
        }
        assert_eq!(
            *warnings.0.lock().unwrap(),
            vec![format!(
                "LoFTEE: the GERP bigWig {} holds no chromosome 22; 22:1500-2000 counts 0 \
                 toward GERP_DIST, as Ensembl VEP's plugin counts it, and later intervals on a \
                 chromosome the bigWig lacks are not reported in this run",
                fixture_gerp_bigwig().display()
            )]
        );
    }

    /// `gerp_dist.pl` opens the bigWig through `Bio::DB::BigWig->new`, which
    /// dies on a path it cannot open, so the row loses its LoF fields under
    /// `Plugin 'LoF' went wrong` and no row is written with a total the file
    /// did not give. A read the library returns as an error after the file
    /// opened fails the row the same way, through `PluginError::Run`, with a
    /// text that names the file and not the interval, so a file failing on
    /// every batch is warned about once.
    #[test]
    fn a_gerp_bigwig_read_failure_fails_the_row() {
        let (p, bigwig, _dir) = plugin_over_a_cut_gerp_bigwig();
        let (tc, v) = gerp_row("21");
        let err = p
            .run(&tc, &v)
            .expect_err("a GERP read failure fails the row");
        let prefix = format!("LoFTEE: reading {}: ", bigwig.display());
        assert!(
            matches!(&err, PluginError::Run(text) if text.starts_with(&prefix) && !text.contains("1500")),
            "{err}"
        );
    }

    /// master's `gerp_dist.pl` reads the tabix TSV through a `tabix` subprocess
    /// with its standard error merged in (`2>&1`) and sums the third field of
    /// every output line, so a read failure's message sums as 0 with no word,
    /// and under master's `+180` cutoff that 0 fires `END_TRUNC` wherever the
    /// 50-base rule fails: a defect, not a count. vep-rs fails the row, as it
    /// fails a bigWig row the file could not give.
    #[test]
    fn a_tabix_gerp_read_failure_fails_the_row() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", "21\t1500\t1.5", "21\t1501\t2.0"],
            IndexSpec::standard(),
        );
        let mut p = LofteePlugin::new();
        p.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        p.init(&[
            format!("gerp_tabix={}", path.display()),
            SPLICE_PREDICTIONS_OFF.to_string(),
        ])
        .unwrap();
        std::fs::remove_file(&path).unwrap();
        let (tc, v) = gerp_row("21");
        let err = p
            .run(&tc, &v)
            .expect_err("a GERP read failure fails the row");
        let prefix = format!("LoFTEE: failed to read {}: ", path.display());
        assert!(
            matches!(&err, PluginError::Run(text) if text.starts_with(&prefix)),
            "{err}"
        );
    }

    /// master's `gerp_dist.pl` reads the tabix TSV through a `tabix`
    /// subprocess, which prints nothing for a region on a sequence its index
    /// does not name, so `get_interval_gerp_tabix` sums 0 toward `GERP_DIST`
    /// without a word, and under master's `+180` cutoff that 0 fires
    /// `END_TRUNC` wherever the 50-base rule fails. vep-rs counts the same 0,
    /// writes the same row and says so once per run, as it does for a bigWig
    /// the chromosome is missing from; a span with no rows on a chromosome the
    /// file holds is a plain 0.
    #[test]
    fn a_chromosome_the_gerp_tabix_tsv_lacks_counts_zero_and_is_warned_once() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec, Warnings};
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", "21\t1500\t1.5", "21\t1501\t2.0"],
            IndexSpec::standard(),
        );
        let mut p = LofteePlugin::new();
        p.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        p.init(&[
            format!("gerp_tabix={}", path.display()),
            SPLICE_PREDICTIONS_OFF.to_string(),
        ])
        .unwrap();
        let warnings = Arc::new(Warnings::default());
        let outputs = tracing::subscriber::with_default(Arc::clone(&warnings), || {
            let (tc, v) = gerp_row("22");
            let first = p.run(&tc, &v).unwrap();
            let (tc, v) = gerp_row("22");
            let second = p.run(&tc, &v).unwrap();
            [first, second]
        });
        for out in &outputs {
            assert_eq!(out.get("LoF").map(String::as_str), Some("LC"), "{out:?}");
            assert_eq!(
                out.get("LoF_filter").map(String::as_str),
                Some("END_TRUNC"),
                "{out:?}"
            );
            assert_eq!(
                out.get("LoF_info").map(String::as_str),
                Some(
                    "GERP_DIST:0,BP_DIST:500,PERCENTILE:0.166666666666667,\
                     DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                )
            );
        }
        assert_eq!(
            *warnings.0.lock().unwrap(),
            vec![format!(
                "LoFTEE: the GERP tabix TSV {} holds no chromosome 22; 22:1500-2000 counts 0 \
                 toward GERP_DIST, as Ensembl VEP's plugin counts it, and later intervals on a \
                 chromosome the tabix TSV lacks are not reported in this run",
                path.display()
            )]
        );

        let reader = GerpReader::open(
            &GerpSource::TabixTsv(path.to_string_lossy().to_string()),
            None,
        )
        .expect("open tabix GERP");
        let mut summer = GerpSummer::of(&reader);
        assert!(matches!(
            summer.sum_interval("22", 1500, 2000),
            Err(GerpReadFailure::NoChromosome)
        ));
        assert_eq!(summer.sum_interval("21", 100, 200).unwrap(), 0.0);
        assert_eq!(summer.sum_interval("21", 1500, 1501).unwrap(), 3.5);
    }

    /// The two branches of `LoF.pm` push the GERP rule's `LoF_info` items in
    /// different orders: master, the plugin of a GRCh37 run, `GERP_DIST`,
    /// `BP_DIST`, then `PERCENTILE` (lines 259-261); the grch38 branch
    /// `PERCENTILE` first (lines 239, 247 and 248). One row over one GERP file
    /// is written in its assembly's order, and in the grch38 branch's when no
    /// assembly is named.
    #[test]
    fn gerp_items_of_lof_info_follow_the_order_of_the_assembly_s_branch() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use vep_core::assembly::Assembly;
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", "21\t1500\t1.5", "21\t1501\t2.0"],
            IndexSpec::standard(),
        );
        let info_under = |assembly: Option<Assembly>| {
            let mut p = LofteePlugin::new();
            p.set_assembly(assembly);
            let mut params = vec![format!("gerp_tabix={}", path.display())];
            if assembly == Some(Assembly::Grch37) {
                params.push(SPLICE_PREDICTIONS_OFF.to_string());
            }
            p.init(&params).unwrap();
            let (tc, v) = gerp_row("21");
            p.run(&tc, &v).unwrap().get("LoF_info").cloned()
        };
        assert_eq!(
            info_under(Some(Assembly::Grch37)).as_deref(),
            Some(
                "GERP_DIST:3.5,BP_DIST:500,PERCENTILE:0.166666666666667,\
                 DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
            )
        );
        for assembly in [Some(Assembly::Grch38), None] {
            assert_eq!(
                info_under(assembly).as_deref(),
                Some(
                    "PERCENTILE:0.166666666666667,GERP_DIST:3.5,BP_DIST:500,\
                     DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                ),
                "{assembly:?}"
            );
        }
    }

    /// The fields of `gerp_row("21")` under `assembly` and `words` beside a
    /// tabix GERP file totalling `gerp` over the row's one interval, with the
    /// warnings `init` logged. A GRCh37 run also passes
    /// [`SPLICE_PREDICTIONS_OFF`], ahead of `words`, without which it is refused.
    fn gerp_row_under(
        assembly: Option<vep_core::assembly::Assembly>,
        gerp: &str,
        words: &[&str],
    ) -> (IndexMap<String, String>, Vec<String>) {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec, Warnings};
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", &format!("21\t1500\t{gerp}")],
            IndexSpec::standard(),
        );
        let mut params = vec![format!("gerp_tabix={}", path.display())];
        if assembly == Some(vep_core::assembly::Assembly::Grch37) {
            params.push(SPLICE_PREDICTIONS_OFF.to_string());
        }
        params.extend(words.iter().map(|word| word.to_string()));
        let warnings = Arc::new(Warnings::default());
        let mut p = LofteePlugin::new();
        p.set_assembly(assembly);
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            p.init(&params)
                .unwrap_or_else(|e| panic!("{assembly:?} {words:?}: {e}"));
        });
        let (tc, v) = gerp_row("21");
        let logged = warnings.0.lock().unwrap().clone();
        (p.run(&tc, &v).unwrap(), logged)
    }

    /// master's `LoF.pm`, the plugin of a GRCh37 run, tests
    /// `use_gerp_end_trunc` nowhere (its one line, 79, is a default that never
    /// runs): its percentile test is commented out (line 253) and every
    /// stop-gained or frameshift row with a CDS end carries the GERP items
    /// (lines 251 and 258-261). The grch38 branch tests the key
    /// (line 240). An off word therefore leaves a GRCh37 row as it is without
    /// the key, under a warning naming it, and selects the percentile rule on
    /// GRCh38 and when no assembly is named.
    #[test]
    fn use_gerp_end_trunc_is_read_on_grch38_alone() {
        use vep_core::assembly::Assembly;
        let gerp_rule = "GERP_DIST:3.5,BP_DIST:500,PERCENTILE:0.166666666666667,\
                         DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL";
        let (without, logged) = gerp_row_under(Some(Assembly::Grch37), "3.5", &[]);
        assert_eq!(without.get("LoF_info").map(String::as_str), Some(gerp_rule));
        assert!(logged.is_empty(), "{logged:?}");
        for word in [
            "use_gerp_end_trunc=0",
            "use_gerp_end_trunc:0",
            "use_gerp_end_trunc=false",
            "use_gerp_end_trunc:FALSE",
        ] {
            let (row, logged) = gerp_row_under(Some(Assembly::Grch37), "3.5", &[word]);
            assert_eq!(row, without, "{word}");
            assert_eq!(
                logged,
                [
                    "LoFTEE: ignoring 'use_gerp_end_trunc' on GRCh37: Ensembl VEP's LoF plugin \
                     for that assembly never reads it and judges END_TRUNC by the GERP rule \
                     at its own cutoff of 180"
                ],
                "{word}"
            );
            for assembly in [Some(Assembly::Grch38), None] {
                let (row, logged) = gerp_row_under(assembly, "3.5", &[word]);
                assert_eq!(
                    row.get("LoF_info").map(String::as_str),
                    Some("PERCENTILE:0.166666666666667"),
                    "{assembly:?} {word}"
                );
                assert!(logged.is_empty(), "{assembly:?} {word}: {logged:?}");
            }
        }
    }

    /// master's `LoF.pm` pushes `END_TRUNC` when the distance from the last
    /// exon is at most 50 and `GERP_DIST` at most 180, a literal (line 268), and
    /// its percentile test, the one reader of `filter_position`, is commented
    /// out (lines 74 and 253). The grch38 branch reads `gerp_end_trunc_cutoff`
    /// (line 254) and `filter_position` (line 241). A GERP total of 250 over a
    /// row inside the last exon is therefore high confidence on GRCh37 whatever
    /// the two keys say, each under its warning, and follows them on GRCh38.
    #[test]
    fn gerp_end_trunc_cutoff_and_filter_position_are_read_on_grch38_alone() {
        use vep_core::assembly::Assembly;
        let call = |row: &IndexMap<String, String>| {
            (
                row.get("LoF").cloned(),
                row.get("LoF_filter").cloned(),
                row.get("LoF_info").cloned(),
            )
        };
        let grch37_row = (
            Some("HC".to_string()),
            None,
            Some(
                "GERP_DIST:250,BP_DIST:500,PERCENTILE:0.166666666666667,\
                 DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                    .to_string(),
            ),
        );
        let (without, logged) = gerp_row_under(Some(Assembly::Grch37), "250", &[]);
        assert_eq!(call(&without), grch37_row);
        assert!(logged.is_empty(), "{logged:?}");
        // master's test is `$gerp_dist <= 180`: a total of 180 is filtered and
        // the next one above it is not.
        for (gerp, filter) in [("180", Some("END_TRUNC")), ("180.001", None)] {
            let (row, _) = gerp_row_under(Some(Assembly::Grch37), gerp, &[]);
            assert_eq!(
                row.get("LoF_filter").map(String::as_str),
                filter,
                "GERP total {gerp}"
            );
        }
        let ignored = |key: &str| {
            format!(
                "LoFTEE: ignoring '{key}' on GRCh37: Ensembl VEP's LoF plugin for that \
                 assembly never reads it and judges END_TRUNC by the GERP rule at its own \
                 cutoff of 180"
            )
        };

        let (row, logged) = gerp_row_under(
            Some(Assembly::Grch37),
            "250",
            &["gerp_end_trunc_cutoff=1000"],
        );
        assert_eq!(call(&row), grch37_row);
        assert_eq!(logged, [ignored("gerp_end_trunc_cutoff")]);
        let (row, logged) = gerp_row_under(
            Some(Assembly::Grch37),
            "250",
            &["use_gerp_end_trunc:0", "filter_position:0.9"],
        );
        assert_eq!(call(&row), grch37_row);
        assert_eq!(
            logged,
            [ignored("use_gerp_end_trunc"), ignored("filter_position")]
        );

        let grch38_gerp_items = "PERCENTILE:0.166666666666667,GERP_DIST:250,BP_DIST:500,\
                                 DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL";
        let (default_cutoff, _) = gerp_row_under(Some(Assembly::Grch38), "250", &[]);
        assert_eq!(
            call(&default_cutoff),
            (
                Some("HC".to_string()),
                None,
                Some(grch38_gerp_items.to_string())
            )
        );
        let (row, logged) = gerp_row_under(
            Some(Assembly::Grch38),
            "250",
            &["gerp_end_trunc_cutoff=1000"],
        );
        assert_eq!(
            call(&row),
            (
                Some("LC".to_string()),
                Some("END_TRUNC".to_string()),
                Some(grch38_gerp_items.to_string())
            )
        );
        assert!(logged.is_empty(), "{logged:?}");
        let (row, logged) = gerp_row_under(
            Some(Assembly::Grch38),
            "250",
            &["use_gerp_end_trunc:0", "filter_position:0.9"],
        );
        assert_eq!(
            call(&row),
            (
                Some("LC".to_string()),
                Some("END_TRUNC".to_string()),
                Some("PERCENTILE:0.166666666666667".to_string())
            )
        );
        assert!(logged.is_empty(), "{logged:?}");
    }

    /// master's `LoF.pm` writes `LoF_info` only while `skip_lof_info` is false by
    /// Perl's string truth (lines 80, 243 and 339), so on GRCh37 any value but
    /// `0` and the empty one drops the field and leaves the call and its flags.
    /// The grch38 branch has no such key: there the word is unknown, logged,
    /// and the field stays.
    #[test]
    fn skip_lof_info_drops_lof_info_on_grch37_alone() {
        use vep_core::assembly::Assembly;
        let grch37_info = "GERP_DIST:3.5,BP_DIST:500,PERCENTILE:0.166666666666667,\
                           DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL";
        let grch38_info = "PERCENTILE:0.166666666666667,GERP_DIST:3.5,BP_DIST:500,\
                           DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL";
        for (word, dropped) in [
            ("skip_lof_info:1", true),
            ("skip_lof_info=1", true),
            ("skip_lof_info:false", true),
            ("skip_lof_info: 0", true),
            ("skip_lof_info:0", false),
            ("skip_lof_info=", false),
        ] {
            let (row, logged) = gerp_row_under(Some(Assembly::Grch37), "3.5", &[word]);
            assert_eq!(row.get("LoF").map(String::as_str), Some("LC"), "{word}");
            assert_eq!(
                row.get("LoF_filter").map(String::as_str),
                Some("END_TRUNC"),
                "{word}"
            );
            assert_eq!(
                row.get("LoF_flags").map(String::as_str),
                Some("SINGLE_EXON"),
                "{word}"
            );
            assert_eq!(
                row.get("LoF_info").map(String::as_str),
                (!dropped).then_some(grch37_info),
                "{word}"
            );
            assert!(logged.is_empty(), "{word}: {logged:?}");

            for assembly in [Some(Assembly::Grch38), None] {
                let (row, logged) = gerp_row_under(assembly, "3.5", &[word]);
                assert_eq!(
                    row.get("LoF_info").map(String::as_str),
                    Some(grch38_info),
                    "{assembly:?} {word}"
                );
                let unknown = if word.contains(':') {
                    format!(
                        "LoFTEE: ignoring unknown parameter 'skip_lof_info' in '{word}': a word \
                         holding ':' is read as key:value at that colon, as Ensembl VEP's LoF \
                         plugin reads it"
                    )
                } else {
                    "LoFTEE: ignoring unknown parameter 'skip_lof_info'".to_string()
                };
                assert_eq!(logged, [unknown], "{assembly:?} {word}");
            }
        }
    }

    /// `LoF.pm` keeps one hash entry per key, its last value (lines 67-72), and
    /// tests `conservation_file` and `run_splice_predictions` once every word
    /// is stored (lines 121, 172 and 274). A value asking for the PhyloCSF
    /// items or the splice-prediction layer is therefore refused only as its
    /// key's last word: an off word after it switches the key off, and an off
    /// word before it does not spare it.
    #[test]
    fn a_refused_key_is_judged_by_its_last_value() {
        let init = |words: &[&str]| {
            let mut params: Vec<String> = words.iter().map(|word| word.to_string()).collect();
            params.push("use_gerp_end_trunc=0".to_string());
            LofteePlugin::new().init(&params)
        };
        for words in [
            [
                "conservation_file:/data/phylocsf.sqlite",
                "conservation_file:false",
            ],
            ["run_splice_predictions:1", "run_splice_predictions:0"],
        ] {
            init(&words).unwrap_or_else(|e| panic!("{words:?}: {e}"));
        }
        for (words, refusal) in [
            (
                [
                    "conservation_file:false",
                    "conservation_file:/data/phylocsf.sqlite",
                ],
                "LoFTEE: conservation_file=/data/phylocsf.sqlite asks for the PhyloCSF items of \
                 LoF_flags and LoF_info, which vep-rs does not compute",
            ),
            (
                ["run_splice_predictions:0", "run_splice_predictions:1"],
                "LoFTEE: run_splice_predictions=1 asks for the splice-prediction layer, which \
                 vep-rs does not implement",
            ),
        ] {
            match init(&words) {
                Err(PluginError::Init(message)) => assert_eq!(message, refusal, "{words:?}"),
                other => panic!("{words:?}: expected PluginError::Init, got {other:?}"),
            }
        }
    }

    /// master's `LoF.pm` names the `tabix` program it runs on `gerp_file` with
    /// `tabix_path` (lines 130, 138 and 141), so on GRCh37 the key is known and,
    /// the GERP TSV being read in-process, changes nothing. The grch38 branch has
    /// no such key: there the word is unknown and logged.
    #[test]
    fn tabix_path_is_a_known_key_on_grch37_alone() {
        use vep_core::assembly::Assembly;
        let (without, _) = gerp_row_under(Some(Assembly::Grch37), "3.5", &[]);
        let (row, logged) = gerp_row_under(
            Some(Assembly::Grch37),
            "3.5",
            &["tabix_path:/usr/local/bin/tabix"],
        );
        assert_eq!(row, without);
        assert!(logged.is_empty(), "{logged:?}");

        let (without, _) = gerp_row_under(Some(Assembly::Grch38), "3.5", &[]);
        for assembly in [Some(Assembly::Grch38), None] {
            let (row, logged) =
                gerp_row_under(assembly, "3.5", &["tabix_path:/usr/local/bin/tabix"]);
            assert_eq!(row, without, "{assembly:?}");
            assert_eq!(
                logged,
                ["LoFTEE: ignoring unknown parameter 'tabix_path' in \
                     'tabix_path:/usr/local/bin/tabix': a word holding ':' is read as key:value \
                     at that colon, as Ensembl VEP's LoF plugin reads it"],
                "{assembly:?}"
            );
        }
    }

    /// Ensembl VEP's `run_plugins` evals each row's plugin call on its own
    /// (`OutputFactory.pm` lines 2225-2244), so a row whose GERP read fails
    /// loses its LoF fields alone. The trait's `run_batch` keeps the rows
    /// before and after it, and the variants after it, and returns the row's
    /// failure for the pipeline's `went wrong` warning.
    #[test]
    fn a_failed_gerp_row_loses_its_lof_fields_alone() {
        let (p, _bigwig, _dir) = plugin_over_a_cut_gerp_bigwig();
        // A splice row reads no GERP: `LoF=HC` with no filter or info.
        let splice_row = || TranscriptConsequence {
            transcript: Some(transcript_of_biotype("protein_coding")),
            consequences: smallvec![Consequence::SpliceDonorVariant],
            loftee_ctx: Some(Box::new(coding_ctx(1))),
            ..TranscriptConsequence::default()
        };
        let (gerp_tc, mut first) = gerp_row("21");
        first.transcript_consequences = vec![splice_row(), gerp_tc, splice_row()];
        let (_, mut second) = gerp_row("21");
        second.transcript_consequences = vec![splice_row()];
        let mut variants = vec![first, second];

        let err = p
            .run_batch(&mut variants)
            .expect_err("the failed row's error is the batch's");
        assert!(matches!(err, PluginError::Run(_)), "{err}");
        let lof_of = |variant: &InputVariant| -> Vec<Option<String>> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.plugin_data.get("LoF").cloned())
                .collect()
        };
        let hc = Some("HC".to_string());
        assert_eq!(lof_of(&variants[0]), vec![hc.clone(), None, hc.clone()]);
        assert!(variants[0].transcript_consequences[1]
            .plugin_data
            .is_empty());
        assert_eq!(lof_of(&variants[1]), vec![hc]);
    }

    #[test]
    fn init_missing_ancestor_fa_errors() {
        let mut p = LofteePlugin::new();
        let err = p
            .init(&[
                "human_ancestor_fa=/nonexistent/human_ancestor.fa".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .expect_err("a missing ancestral-allele FASTA must fail at init")
            .to_string();
        assert!(
            err.contains("failed to open human_ancestor_fa"),
            "unexpected error: {err}"
        );
    }

    /// `human_ancestor_fa` opens a BGZF FASTA through its `.fai` and `.gzi`, the
    /// form the upstream bundle ships; without the `.gzi` the open fails and the
    /// error names both accepted forms.
    #[test]
    fn init_opens_a_bgzf_ancestor_fasta() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh37-hgvs");
        let bgzf = root.join("reference.fa.gz");
        let mut p = LofteePlugin::new();
        p.init(&[
            format!("human_ancestor_fa={}", bgzf.display()),
            "use_gerp_end_trunc=0".to_string(),
        ])
        .unwrap();
        assert!(p.ancestor_fasta.is_some());

        let dir = tempfile::tempdir().unwrap();
        for name in ["reference.fa.gz", "reference.fa.gz.fai"] {
            std::fs::copy(root.join(name), dir.path().join(name)).unwrap();
        }
        let mut without_gzi = LofteePlugin::new();
        let err = without_gzi
            .init(&[
                format!(
                    "human_ancestor_fa={}",
                    dir.path().join("reference.fa.gz").display()
                ),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("a plain .fa with its .fai, or a bgzipped .fa.gz with its .fai and .gzi"),
            "{err}"
        );
    }

    #[test]
    fn init_accepts_colon_syntax() {
        let mut p = LofteePlugin::new();
        p.init(&[
            "min_intron_size:25".to_string(),
            "use_gerp_end_trunc:0".to_string(),
        ])
        .unwrap();
        assert_eq!(p.min_intron_size, 25.0);
    }

    /// A word holding a colon is LoFTEE's `key:value`, read whatever `=` it
    /// holds (`LoF.pm` lines 67-72 split on `:` alone), so a `loftee_path` whose
    /// directory name holds `=` derives the GERP bigWig from that directory; a
    /// key:value word whose key holds `=` is the unknown key Ensembl VEP's
    /// plugin stores and never reads. `LoF.pm` line 69 keeps a word only as
    /// exactly two fields, so a value holding a second colon (`key:value:x`,
    /// `key::value`) is dropped with the whole word and the key keeps its
    /// default. Trailing colons are dropped as Perl's `split` drops them, so
    /// `key:value:` is read while `key:` and `key=value:` are not. A word
    /// without a colon is vep-rs's `key=value`, read at its first `=`.
    #[test]
    fn parameters_split_on_the_colon_alone() {
        for (word, split) in [
            ("key:value", Some(("key", "value"))),
            ("key:value:", Some(("key", "value"))),
            ("key:value::", Some(("key", "value"))),
            (":value", Some(("", "value"))),
            ("key=value:x", Some(("key=value", "x"))),
            ("key:value:x", None),
            ("key::value", None),
            ("::value", None),
            ("key:", None),
            (":", None),
            ("gerp_bigwig:https://example.invalid/gerp.bw", None),
            ("key=value", Some(("key", "value"))),
            ("key=a=b", Some(("key", "a=b"))),
            ("key", None),
        ] {
            assert_eq!(split_parameter(word), split, "{word}");
        }

        let dir = tempfile::tempdir().unwrap();
        let loftee_path = dir.path().join("loftee=a46b502a");
        std::fs::create_dir(&loftee_path).unwrap();
        let bigwig = loftee_path.join(GRCH38_DEFAULT_GERP_BIGWIG);
        std::fs::copy(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(
                "../../tests/golden/116/GRCh38-plugin-loftee/plugin_data/\
                 gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            ),
            &bigwig,
        )
        .unwrap();

        let loftee_path_value = format!("{}/", loftee_path.display());
        let mut derived = LofteePlugin::new();
        derived
            .init(&[format!("loftee_path:{loftee_path_value}")])
            .unwrap();
        assert_eq!(
            derived.loftee_path.as_deref(),
            Some(loftee_path_value.as_str())
        );
        assert!(
            matches!(&derived.gerp_source, Some(GerpSource::BigWig(p)) if *p == bigwig.display().to_string()),
            "{:?}",
            derived.gerp_source
        );

        let mut second_colon = LofteePlugin::new();
        second_colon
            .init(&[
                "loftee_path:/data/a:b/".to_string(),
                "gerp_bigwig:/data/a:b/gerp.bw".to_string(),
                "min_intron_size:20:5".to_string(),
                "human_ancestor_fa::/data/anc.fa.gz".to_string(),
                "use_gerp_end_trunc:0".to_string(),
            ])
            .unwrap();
        assert_eq!(second_colon.loftee_path, None);
        assert!(second_colon.gerp_source.is_none());
        assert_eq!(second_colon.min_intron_size, DEFAULT_MIN_INTRON_SIZE);
        assert_eq!(second_colon.human_ancestor_fa, None);
        assert!(!second_colon.use_gerp_end_trunc);

        let mut equals_in_key = LofteePlugin::new();
        equals_in_key
            .init(&[
                "gerp_bigwig=/data/x:y.bw".to_string(),
                "use_gerp_end_trunc:0".to_string(),
            ])
            .unwrap();
        assert!(equals_in_key.gerp_source.is_none());

        let mut trailing_colons = LofteePlugin::new();
        trailing_colons
            .init(&[
                "min_intron_size:20:".to_string(),
                "filter_position:".to_string(),
                "gerp_bigwig=/data/gerp.bw:".to_string(),
                "use_gerp_end_trunc:0::".to_string(),
            ])
            .unwrap();
        assert_eq!(trailing_colons.min_intron_size, 20.0);
        assert_eq!(trailing_colons.filter_position, DEFAULT_FILTER_POSITION);
        assert!(trailing_colons.gerp_source.is_none());
        assert!(!trailing_colons.use_gerp_end_trunc);

        let mut plain = LofteePlugin::new();
        plain
            .init(&[
                "min_intron_size:25".to_string(),
                "filter_position=0.1".to_string(),
                "loftee_path=/data/loftee=a46b502a/".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        assert_eq!(plain.min_intron_size, 25.0);
        assert_eq!(plain.filter_position, 0.1);
        assert_eq!(plain.loftee_path.as_deref(), Some("/data/loftee=a46b502a/"));
    }

    /// A `key:value` word whose key Ensembl VEP's LoF plugin stores and never
    /// reads is dropped with a warning naming the key and the whole word, as
    /// the key alone does not show the word it came from: `gerp_bigwig=/data/x:y.bw`
    /// splits at its colon into the key `gerp_bigwig=/data/x`, and a bare URL
    /// splits at the colon of its scheme into the scheme alone.
    #[test]
    fn an_unknown_colon_key_is_warned_with_its_whole_word() {
        use crate::test_fixtures::Warnings;
        let warnings = Arc::new(Warnings::default());
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            LofteePlugin::new()
                .init(&[
                    "gerp_bigwig=/data/x:y.bw".to_string(),
                    "https://example.invalid/gerp.bw".to_string(),
                    "use_gerp_end_trunc=0".to_string(),
                ])
                .unwrap();
        });
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                "LoFTEE: ignoring unknown parameter 'gerp_bigwig=/data/x' in \
                 'gerp_bigwig=/data/x:y.bw': a word holding ':' is read as key:value at that \
                 colon, as Ensembl VEP's LoF plugin reads it",
                "LoFTEE: ignoring unknown parameter 'https' in 'https://example.invalid/gerp.bw': \
                 a word holding ':' is read as key:value at that colon, as Ensembl VEP's LoF \
                 plugin reads it",
            ]
        );
    }

    /// `run_splice_predictions` is read by Perl's string truth (LoF.pm:172), so
    /// `1` and `true` ask for the splice-prediction layer, which is refused at
    /// initialization, while `0` leaves it off in either grammar, and so does
    /// an empty value in a run naming no assembly, as this one.
    #[test]
    fn run_splice_predictions_true_is_refused() {
        for value in ["1", "true"] {
            let mut p = LofteePlugin::new();
            let err = p
                .init(&[
                    format!("run_splice_predictions={value}"),
                    "use_gerp_end_trunc=0".to_string(),
                ])
                .expect_err("a requested splice-prediction layer is refused");
            match err {
                PluginError::Init(message) => assert_eq!(
                    message,
                    format!("LoFTEE: run_splice_predictions={value} asks for the splice-prediction layer, which vep-rs does not implement")
                ),
                other => panic!("expected PluginError::Init, got {other:?}"),
            }
        }
        for param in [
            "run_splice_predictions=0",
            "run_splice_predictions:0",
            "run_splice_predictions=",
        ] {
            LofteePlugin::new()
                .init(&[param.to_string(), "use_gerp_end_trunc=0".to_string()])
                .unwrap();
        }
    }

    /// `LoF.pm` keeps a field as `split /:/` yields it (lines 68-70), so a
    /// value padded with a blank is another string: ` 0` is true at the
    /// `run_splice_predictions` truth test (line 172), so the word asks for
    /// the splice-prediction layer, as it does there.
    #[test]
    fn run_splice_predictions_with_a_blank_before_its_zero_is_true_as_in_perl() {
        for param in ["run_splice_predictions: 0", "run_splice_predictions:0 "] {
            let err = LofteePlugin::new()
                .init(&[param.to_string(), "use_gerp_end_trunc=0".to_string()])
                .expect_err(param);
            let value = param.split_once(':').unwrap().1;
            match err {
                PluginError::Init(message) => assert_eq!(
                    message,
                    format!(
                        "LoFTEE: run_splice_predictions={value} asks for the splice-prediction \
                         layer, which vep-rs does not implement; the padded value '{value}' is \
                         true in Ensembl VEP's plugin"
                    ),
                    "{param}"
                ),
                other => panic!("{param}: expected PluginError::Init, got {other:?}"),
            }
        }
    }

    /// master's `LoF.pm`, the plugin of a GRCh37 run, defaults
    /// `run_splice_predictions` to 1 (line 98) and tests it by Perl's string
    /// truth (line 186), so the splice-prediction layer runs there unless the
    /// key's last stored value is `0`; the grch38 branch defaults the key to 0
    /// (its line 98). A GRCh37 run whose key is absent, empty, or left at its
    /// default by a one-field word `split /:/` drops is therefore refused at
    /// initialization, naming the word to pass. The same words start a run on
    /// GRCh38 and when no assembly is named, and a true value is refused on
    /// every assembly.
    #[test]
    fn a_grch37_run_that_leaves_splice_predictions_on_is_refused() {
        use vep_core::assembly::Assembly;
        let dir = tempfile::tempdir().unwrap();
        let gerp = format!("gerp_tabix={}", fixture_gerp_tabix(dir.path()).display());
        let init = |assembly: Option<Assembly>, words: &[&str]| {
            let mut params = vec![gerp.clone()];
            params.extend(words.iter().map(|word| word.to_string()));
            let mut p = LofteePlugin::new();
            p.set_assembly(assembly);
            p.init(&params)
        };
        let left_on: [&[&str]; 4] = [
            &[],
            &["run_splice_predictions="],
            &["run_splice_predictions:"],
            &[SPLICE_PREDICTIONS_OFF, "run_splice_predictions="],
        ];
        let refusal = "LoFTEE: this GRCh37 run asks for the splice-prediction layer, which vep-rs \
                       does not implement: Ensembl VEP's LoF plugin runs it on that assembly \
                       unless run_splice_predictions is 0; pass run_splice_predictions:0";
        for words in left_on {
            match init(Some(Assembly::Grch37), words) {
                Err(PluginError::Init(message)) => assert_eq!(message, refusal, "{words:?}"),
                other => panic!("{words:?}: expected PluginError::Init, got {other:?}"),
            }
            for assembly in [Some(Assembly::Grch38), None] {
                init(assembly, words).unwrap_or_else(|e| panic!("{assembly:?} {words:?}: {e}"));
            }
        }
        // The key is judged before the GERP source opens, so a run lacking the
        // key and naming a file that is not there is refused for the key.
        let mut no_gerp_file = LofteePlugin::new();
        no_gerp_file.set_assembly(Some(Assembly::Grch37));
        match no_gerp_file.init(&["gerp_file:/nonexistent/GERP_scores.txt.gz".to_string()]) {
            Err(PluginError::Init(message)) => assert_eq!(message, refusal),
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
        let switched_off: [&[&str]; 4] = [
            &[SPLICE_PREDICTIONS_OFF],
            &["run_splice_predictions=0"],
            &["run_splice_predictions:1", SPLICE_PREDICTIONS_OFF],
            &[SPLICE_PREDICTIONS_OFF, "run_splice_predictions:"],
        ];
        for words in switched_off {
            init(Some(Assembly::Grch37), words).unwrap_or_else(|e| panic!("{words:?}: {e}"));
        }
        // A value true by Perl's string truth, a zero beside a blank among them
        // (`LoF.pm` lines 68-70 keep the blank), asks for the layer on every
        // assembly.
        for assembly in [Some(Assembly::Grch37), Some(Assembly::Grch38), None] {
            for (word, refusal) in [
                (
                    "run_splice_predictions:1",
                    "LoFTEE: run_splice_predictions=1 asks for the splice-prediction layer, \
                     which vep-rs does not implement",
                ),
                (
                    "run_splice_predictions: 0",
                    "LoFTEE: run_splice_predictions= 0 asks for the splice-prediction layer, \
                     which vep-rs does not implement; the padded value ' 0' is true in Ensembl \
                     VEP's plugin",
                ),
            ] {
                match init(assembly, &[SPLICE_PREDICTIONS_OFF, word]) {
                    Err(PluginError::Init(message)) => {
                        assert_eq!(message, refusal, "{assembly:?} {word}")
                    }
                    other => {
                        panic!("{assembly:?} {word}: expected PluginError::Init, got {other:?}")
                    }
                }
            }
        }
    }

    /// The word that switches the splice-prediction layer off moves no field:
    /// a GRCh37 run, which must pass it, writes master's row, and a run on
    /// GRCh38 or naming no assembly writes one row with the word, with an
    /// empty value and without the key, the grch38 branch's default being 0
    /// (`LoF.pm` line 98).
    #[test]
    fn run_splice_predictions_off_moves_no_field() {
        use vep_core::assembly::Assembly;
        fn fields(row: &IndexMap<String, String>) -> Vec<(&str, &str)> {
            row.iter()
                .map(|(field, value)| (field.as_str(), value.as_str()))
                .collect()
        }
        let (grch37, logged) = gerp_row_under(Some(Assembly::Grch37), "3.5", &[]);
        assert_eq!(
            fields(&grch37),
            [
                ("LoF", "LC"),
                ("LoF_filter", "END_TRUNC"),
                ("LoF_flags", "SINGLE_EXON"),
                (
                    "LoF_info",
                    "GERP_DIST:3.5,BP_DIST:500,PERCENTILE:0.166666666666667,\
                     DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                ),
            ]
        );
        assert!(logged.is_empty(), "{logged:?}");
        for assembly in [Some(Assembly::Grch38), None] {
            let (without, logged) = gerp_row_under(assembly, "3.5", &[]);
            assert_eq!(
                fields(&without),
                [
                    ("LoF", "HC"),
                    ("LoF_flags", "SINGLE_EXON"),
                    (
                        "LoF_info",
                        "PERCENTILE:0.166666666666667,GERP_DIST:3.5,BP_DIST:500,\
                         DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                    ),
                ],
                "{assembly:?}"
            );
            assert!(logged.is_empty(), "{assembly:?}: {logged:?}");
            for word in [
                SPLICE_PREDICTIONS_OFF,
                "run_splice_predictions=0",
                "run_splice_predictions=",
            ] {
                let (row, logged) = gerp_row_under(assembly, "3.5", &[word]);
                assert_eq!(row, without, "{assembly:?} {word}");
                assert!(logged.is_empty(), "{assembly:?} {word}: {logged:?}");
            }
        }
    }

    /// A key padded with a blank is another hash key in `LoF.pm`, stored and
    /// never read, so the plugin keeps the default and vep-rs logs the word.
    #[test]
    fn a_key_with_a_blank_is_another_key_as_in_perl() {
        use crate::test_fixtures::Warnings;
        let warnings = Arc::new(Warnings::default());
        let mut p = LofteePlugin::new();
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            p.init(&[
                " run_splice_predictions:1".to_string(),
                " min_intron_size:20".to_string(),
                "min_intron_size :20".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        });
        assert_eq!(p.min_intron_size, DEFAULT_MIN_INTRON_SIZE);
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [
                "LoFTEE: ignoring unknown parameter ' run_splice_predictions' in \
                 ' run_splice_predictions:1': a word holding ':' is read as key:value at that \
                 colon, as Ensembl VEP's LoF plugin reads it",
                "LoFTEE: ignoring unknown parameter ' min_intron_size' in ' min_intron_size:20': \
                 a word holding ':' is read as key:value at that colon, as Ensembl VEP's LoF \
                 plugin reads it",
                "LoFTEE: ignoring unknown parameter 'min_intron_size ' in 'min_intron_size :20': \
                 a word holding ':' is read as key:value at that colon, as Ensembl VEP's LoF \
                 plugin reads it",
            ]
        );
    }

    /// The numeric parameters are read where `LoF.pm` compares them as numbers
    /// (`min_intron_size` line 297, `filter_position` line 241,
    /// `gerp_end_trunc_cutoff` line 254), and Perl's numeric reading skips the
    /// blanks around a number (Perl's `isSPACE` set, a vertical tab and a form
    /// feed included) and takes the leading number of the text, so ` 20` and
    /// `20 ` are 20 under either grammar, ` 20abc` is 20 and `abc ` is 0.
    #[test]
    fn a_numeric_value_reads_through_its_blanks_as_perl_numifies_it() {
        for (params, filter_position) in [
            (
                [
                    "min_intron_size: 20",
                    "filter_position: 0.1",
                    "gerp_end_trunc_cutoff: -60",
                ],
                0.1,
            ),
            (
                [
                    "min_intron_size=20 ",
                    "filter_position=\t0.1\t",
                    "gerp_end_trunc_cutoff=-60 ",
                ],
                0.1,
            ),
            (
                [
                    "min_intron_size: 20abc",
                    "filter_position:abc ",
                    "gerp_end_trunc_cutoff:\t-60x",
                ],
                0.0,
            ),
            (
                [
                    "min_intron_size:\x0b20",
                    "filter_position:\x0c0.1\x0b",
                    "gerp_end_trunc_cutoff:-60\r",
                ],
                0.1,
            ),
        ] {
            let mut p = LofteePlugin::new();
            let mut words: Vec<String> = params.iter().map(|w| w.to_string()).collect();
            words.push("use_gerp_end_trunc=0".to_string());
            p.init(&words).unwrap_or_else(|e| panic!("{params:?}: {e}"));
            assert_eq!(p.min_intron_size, 20.0, "{params:?}");
            assert_eq!(p.filter_position, filter_position, "{params:?}");
            assert_eq!(p.gerp_end_trunc_cutoff, -60.0, "{params:?}");
            assert!(p.gerp_cutoff_explicit, "{params:?}");
        }
    }

    /// `use_gerp_end_trunc` is tested for Perl truth (`LoF.pm` line 240) and
    /// `check_complete_cds` compared with `'true'` after `lc` (line 268), so a
    /// blank-padded ` 0` keeps the GERP rule on and a padded ` true` leaves
    /// `INCOMPLETE_CDS` off, as both do there.
    #[test]
    fn a_truth_or_word_tested_value_keeps_its_blank_as_in_perl() {
        let mut gerp_on = LofteePlugin::new();
        let result = gerp_on.init(&[
            "use_gerp_end_trunc: 0".to_string(),
            "check_complete_cds: true".to_string(),
            "gerp_bigwig=/nonexistent/gerp.bw".to_string(),
        ]);
        assert!(gerp_on.use_gerp_end_trunc);
        assert!(!gerp_on.check_complete_cds);
        assert!(
            result.is_err(),
            "the GERP rule stays on and opens its source"
        );
    }

    /// `master`'s `LoF.pm` runs `tabix` on `gerp_file` as a command Perl splits
    /// into words itself (lines 138 and 141), which drops the blanks around the
    /// path, so a padded tabix path opens; the `grch38` branch opens its bigWig
    /// by the name given (`gerp_dist.pl` line 72), so a padded bigWig path is a
    /// missing file.
    #[test]
    fn a_padded_gerp_path_opens_where_perl_runs_a_command_on_it() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let dir = tempfile::tempdir().unwrap();
        let tabix = write_tabix_fixture(
            dir.path(),
            "gerp_scores.txt.gz",
            &["#chrom\tpos\tGERP", "21\t1500\t1.5"],
            IndexSpec::standard(),
        );
        let word = format!(" {}\t", tabix.display());
        let mut padded = LofteePlugin::new();
        padded.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        padded
            .init(&[
                format!("gerp_file:{word}"),
                SPLICE_PREDICTIONS_OFF.to_string(),
            ])
            .unwrap();
        assert_eq!(
            padded.gerp_source.as_ref().map(GerpSource::path),
            Some(word.as_str())
        );
        assert_eq!(
            padded.gerp_source.as_ref().map(GerpSource::opened_path),
            Some(tabix.to_str().unwrap())
        );
        assert!(padded.gerp_reader.is_some());

        let mut missing = LofteePlugin::new();
        missing.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        let err = missing
            .init(&[
                "gerp_file: false".to_string(),
                SPLICE_PREDICTIONS_OFF.to_string(),
            ])
            .expect_err("a padded ` false` is a path, not the off word")
            .to_string();
        assert!(
            err.contains("GERP data file not found: ' false'"),
            "unexpected error: {err}"
        );

        let bigwig = fixture_gerp_bigwig();
        let mut as_given = LofteePlugin::new();
        as_given.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        let err = as_given
            .init(&[format!("gerp_bigwig: {}", bigwig.display())])
            .expect_err("a padded bigWig path is the file named")
            .to_string();
        assert!(
            err.contains(&format!("' {}'", bigwig.display())),
            "unexpected error: {err}"
        );
    }

    /// `LoF.pm` line 312 tests the `human_ancestor_fa` value whole, so a padded
    /// ` false` is a path there, then line 539 runs `samtools faidx` on it as a
    /// command Perl splits into words itself, which drops the blanks around the
    /// path: ` false` is the missing file `false`, and a padded real path opens.
    #[test]
    fn human_ancestor_fa_opens_the_path_perl_hands_samtools() {
        let err = LofteePlugin::new()
            .init(&[
                "human_ancestor_fa: false".to_string(),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .expect_err("a padded ` false` is a path, not the off word")
            .to_string();
        assert!(
            err.contains("failed to open human_ancestor_fa ' false'"),
            "unexpected error: {err}"
        );

        let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-plugin-loftee/plugin_data/human_ancestor.fa.gz");
        let fasta = fasta.to_str().unwrap();
        let word = format!(" {fasta}\t");
        let mut padded = LofteePlugin::new();
        padded
            .init(&[
                format!("human_ancestor_fa:{word}"),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
        assert_eq!(padded.human_ancestor_fa.as_deref(), Some(word.as_str()));
        assert!(padded.ancestor_fasta.is_some());
    }

    /// `LoF.pm` does not open `loftee_path` itself: the value is joined to the
    /// file names under it as given (lines 86-116 and 130), so it keeps its
    /// blanks, and the GERP bigWig derived from a padded `loftee_path` carries
    /// them.
    #[test]
    fn a_padded_loftee_path_stays_as_given() {
        let mut p = LofteePlugin::new();
        p.init(&[
            "loftee_path: /vep/loftee/".to_string(),
            "use_gerp_end_trunc=0".to_string(),
        ])
        .unwrap();
        assert_eq!(p.loftee_path.as_deref(), Some(" /vep/loftee/"));
        assert_eq!(
            p.default_gerp_source().path(),
            format!(" /vep/loftee/{GRCH38_DEFAULT_GERP_BIGWIG}")
        );
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

    /// An empty transcript of `biotype`.
    fn transcript_of_biotype(biotype: &str) -> Arc<Transcript> {
        Arc::new(Transcript {
            biotype: biotype.into(),
            ..Transcript::default()
        })
    }

    fn lof_variant(
        consequences: ConsequenceList,
        ctx: LofteeContext,
    ) -> (TranscriptConsequence, InputVariant) {
        let tc = TranscriptConsequence {
            transcript: Some(transcript_of_biotype("protein_coding")),
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
            transcript: Some(transcript_of_biotype("lncRNA")),
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
            transcript: Some(transcript_of_biotype("protein_coding")),
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

    /// Ensembl VEP 116.2 never hands LoFTEE a structural variant: `LoF.pm`
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2220). A `<DEL>` whose transcript
    /// row carries a stop-gained call gets no LoF field through the batch
    /// dispatch, while the same row on a sequence variant is called.
    #[test]
    fn a_structural_variant_row_gets_no_lof_call() {
        let p = LofteePlugin::new();
        assert!(!p.annotates_structural_variants());
        let row = || TranscriptConsequence {
            transcript: Some(transcript_of_biotype("protein_coding")),
            consequences: smallvec![Consequence::StopGained],
            exon: Some("1/1".to_string()),
            cds_position: Some("150".to_string()),
            loftee_ctx: Some(Box::new(coding_ctx(1))),
            ..TranscriptConsequence::default()
        };
        let mut deletion =
            InputVariant::new("21".into(), 1500, 1600, b"N".to_vec(), b"<DEL>".to_vec());
        deletion.is_structural = true;
        deletion.transcript_consequences.push(row());
        let mut sequence = InputVariant::new("21".into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        sequence.transcript_consequences.push(row());
        let mut variants = vec![deletion, sequence];
        p.run_batch(&mut variants).unwrap();
        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
        assert_eq!(
            variants[1].transcript_consequences[0]
                .plugin_data
                .get("LoF")
                .map(String::as_str),
            Some("HC")
        );
    }

    /// LoF.pm runs its exon checks only when the consequence carries an exon
    /// number, so a frameshift without one (an intronic combination) keeps its
    /// high-confidence call instead of an EXON_INTRON_UNDEF filter.
    #[test]
    fn no_exon_number_keeps_the_call_high_confidence() {
        let mut p = LofteePlugin::new();
        p.init(&["use_gerp_end_trunc=false".to_string()]).unwrap();
        let mut ctx = coding_ctx(1);
        ctx.introns = vec![vep_core::transcript::Intron {
            start: 1400,
            end: 1499,
            rank: 2,
        }];
        let tc = TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                biotype: "protein_coding".into(),
                ..Transcript::default()
            })),
            consequences: smallvec![Consequence::FrameshiftVariant, Consequence::IntronVariant],
            exon: None,
            intron: Some("2/4".to_string()),
            cds_position: Some("150".to_string()),
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        };
        let v = InputVariant::new("21".into(), 1500, 1500, b"C".to_vec(), b"T".to_vec());
        let out = p.run(&tc, &v).unwrap();
        assert_eq!(out.get("LoF").map(String::as_str), Some("HC"));
        assert!(out.get("LoF_filter").is_none(), "{out:?}");
        assert_eq!(
            out.get("LoF_info").map(String::as_str),
            Some("PERCENTILE:0.166666666666667,INTRON_SIZE:100")
        );
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
            transcript: Some(transcript_of_biotype("protein_coding")),
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

    /// A splice-donor row over intron 2 of `size` bases, the row of
    /// `small_intron_splice_is_lc` with the intron's size chosen.
    fn splice_donor_over_an_intron_of(size: u64) -> (TranscriptConsequence, InputVariant) {
        let mut ctx = coding_ctx(1);
        ctx.introns = vec![vep_core::transcript::Intron {
            start: 1490,
            end: 1490 + size - 1,
            rank: 2,
        }];
        let tc = TranscriptConsequence {
            transcript: Some(transcript_of_biotype("protein_coding")),
            consequences: smallvec![Consequence::SpliceDonorVariant],
            intron: Some("2/4".to_string()),
            loftee_ctx: Some(Box::new(ctx)),
            ..TranscriptConsequence::default()
        };
        let v = InputVariant::new("21".into(), 1491, 1491, b"C".to_vec(), b"T".to_vec());
        (tc, v)
    }

    /// `LoF.pm` reads `min_intron_size` in `$intron_size < $self->{min_intron_size}`
    /// alone (line 297), which numifies the text: `15.0` is 15, `15.5` keeps its
    /// fraction, and `abc`, `-1` and `0.05x` are 0, -1 and 0.05, below which no
    /// intron's size falls, so `SMALL_INTRON` is never flagged.
    #[test]
    fn min_intron_size_is_read_as_its_numeric_comparison_reads_it() {
        let filter_of = |value: &str, intron_size: u64| {
            let mut p = LofteePlugin::new();
            p.init(&[
                format!("min_intron_size:{value}"),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
            let (tc, v) = splice_donor_over_an_intron_of(intron_size);
            let out = p.run(&tc, &v).unwrap();
            assert_eq!(
                out.get("LoF_info").map(String::as_str),
                Some(format!("INTRON_SIZE:{intron_size}").as_str()),
                "{value} over {intron_size}"
            );
            let filter = out.get("LoF_filter").cloned();
            let confidence = if filter.is_some() { "LC" } else { "HC" };
            assert_eq!(
                out.get("LoF").map(String::as_str),
                Some(confidence),
                "{value} over {intron_size}"
            );
            filter
        };
        let small_intron = Some("SMALL_INTRON".to_string());
        assert_eq!(filter_of("15.0", 14), small_intron);
        assert_eq!(filter_of("15.0", 15), None);
        assert_eq!(filter_of("15.5", 15), small_intron);
        assert_eq!(filter_of("15.5", 16), None);
        for value in ["abc", "-1", "0.05x"] {
            assert_eq!(filter_of(value, 1), None, "{value}");
            assert_eq!(filter_of(value, 6), None, "{value}");
        }
    }

    /// A stop-gained row in exon 2 of 5 at `cds_position` of the 900-base CDS
    /// of `coding_ctx`.
    fn stop_gained_at_cds_position(cds_position: &str) -> (TranscriptConsequence, InputVariant) {
        let (mut tc, v) = lof_variant(smallvec![Consequence::StopGained], coding_ctx(1));
        tc.cds_position = Some(cds_position.to_string());
        (tc, v)
    }

    /// `LoF.pm` reads `filter_position` in `$lof_percentile >= 1 -
    /// $self->{filter_position}` alone (line 241), which numifies the text:
    /// `0.05x` is 0.05, `abc` is 0, under which `END_TRUNC` needs the whole
    /// CDS, percentile 1, and `-1` is -1, under which no percentile reaches 2.
    #[test]
    fn filter_position_is_read_as_its_arithmetic_reads_it() {
        let filter_of = |value: &str, cds_position: &str| {
            let mut p = LofteePlugin::new();
            p.init(&[
                format!("filter_position:{value}"),
                "use_gerp_end_trunc=0".to_string(),
            ])
            .unwrap();
            let (tc, v) = stop_gained_at_cds_position(cds_position);
            p.run(&tc, &v).unwrap().get("LoF_filter").cloned()
        };
        let end_trunc = Some("END_TRUNC".to_string());
        assert_eq!(filter_of("abc", "900"), end_trunc);
        assert_eq!(filter_of("abc", "899"), None);
        assert_eq!(filter_of("0.05x", "899"), end_trunc);
        assert_eq!(filter_of("0.05x", "854"), None);
        assert_eq!(filter_of("-1", "900"), None);
    }

    /// `LoF.pm` reads `gerp_end_trunc_cutoff` in `$gerp_dist <=
    /// $self->{gerp_end_trunc_cutoff}` alone (line 254), which numifies the
    /// text: `abc` is 0, `0.05x` is 0.05 and `-1` is -1. A chromosome the bigWig
    /// lacks totals a `GERP_DIST` of 0, so a row within 50 bases of its last
    /// exon is `END_TRUNC` under 0 or 0.05 and not under -1 or -58.0. A given
    /// value stands over the assembly's default, as a defined hash entry does
    /// there (line 78).
    #[test]
    fn gerp_end_trunc_cutoff_is_read_as_its_numeric_comparison_reads_it() {
        let filter_of = |value: &str| {
            let mut p = LofteePlugin::new();
            p.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
            p.init(&[
                format!("gerp_end_trunc_cutoff:{value}"),
                format!("gerp_bigwig={}", fixture_gerp_bigwig().display()),
            ])
            .unwrap();
            assert!(p.gerp_cutoff_explicit, "{value}");
            let (tc, v) = gerp_row("22");
            let out = p.run(&tc, &v).unwrap();
            assert_eq!(
                out.get("LoF_info").map(String::as_str),
                Some(
                    "PERCENTILE:0.166666666666667,GERP_DIST:0,BP_DIST:500,\
                     DIST_FROM_LAST_EXON:-500,50_BP_RULE:FAIL"
                ),
                "{value}"
            );
            out.get("LoF_filter").cloned()
        };
        let end_trunc = Some("END_TRUNC".to_string());
        assert_eq!(filter_of("abc"), end_trunc);
        assert_eq!(filter_of("0.05x"), end_trunc);
        assert_eq!(filter_of("-1"), None);
        assert_eq!(filter_of("-58.0"), None);
    }

    /// `GC_TO_GT_DONOR` reads `feature_seq`, the alleles on the transcript's
    /// strand. Against a forward-strand transcript whose intron 2 opens `GC` at
    /// 1490-1491, the donor SNV `C/T` at 1491 strengthens the donor, and the same
    /// change given on the reverse strand (`G/A -1`) is the same `C/T` on the
    /// transcript, while a `C/T -1` reads `G/A` and does not. Against a
    /// reverse-strand transcript whose intron 2 opens `GC` at 1581-1580 (the
    /// genome's `GC`), the donor SNV at 1580 is `G/A` on the genome: given on the
    /// reverse strand it is written `C/T`, read back as `C/T` on the transcript,
    /// and strengthens the donor, while `G/A -1` does not.
    #[test]
    fn gc_to_gt_donor_reads_a_reverse_strand_record_on_the_forward_strand() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let mut sequence = vec![b'A'; 1600];
        sequence[1489] = b'G';
        sequence[1490] = b'C';
        sequence[1579] = b'G';
        sequence[1580] = b'C';
        let mut contents = b">21\n".to_vec();
        contents.extend_from_slice(&sequence);
        contents.push(b'\n');
        std::fs::write(&fa, contents).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t1600\t4\t1600\t1601\n").unwrap();
        let mut p = LofteePlugin::new();
        p.set_reference_fasta(Some(Arc::new(
            vep_fasta::IndexedFasta::from_path(&fa).unwrap(),
        )));
        let transcript = |strand: i8, intron_end: u64| {
            let mut ctx = coding_ctx(strand);
            ctx.introns = vec![vep_core::transcript::Intron {
                start: 1490,
                end: intron_end,
                rank: 2,
            }];
            ctx
        };
        let snv = |position: u64, reference: &[u8], alt: &[u8], strand: Strand| {
            let mut v = InputVariant::new(
                "21".into(),
                position,
                position,
                reference.to_vec(),
                alt.to_vec(),
            );
            v.strand = strand;
            v
        };
        let donor = |ctx: &LofteeContext, v: &InputVariant| p.gc_to_gt_donor(ctx, v, "21", 2);
        let forward = transcript(1, 1580);
        assert_eq!(
            donor(&forward, &snv(1491, b"C", b"T", Strand::Forward)),
            Some(true)
        );
        assert_eq!(
            donor(&forward, &snv(1491, b"G", b"A", Strand::Reverse)),
            Some(true)
        );
        assert_eq!(
            donor(&forward, &snv(1491, b"C", b"T", Strand::Reverse)),
            Some(false)
        );
        let reverse = transcript(-1, 1581);
        assert_eq!(
            donor(&reverse, &snv(1580, b"G", b"A", Strand::Forward)),
            Some(true)
        );
        assert_eq!(
            donor(&reverse, &snv(1580, b"C", b"T", Strand::Reverse)),
            Some(true)
        );
        assert_eq!(
            donor(&reverse, &snv(1580, b"G", b"A", Strand::Reverse)),
            Some(false)
        );
    }

    /// `nagnag_site` for `variant` against a transcript on `strand`, with
    /// `context` the genome's bases 1496 to 1504 and `A` everywhere else.
    fn nagnag_probe_of(context: &[u8; 9], strand: i8, variant: &InputVariant) -> bool {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let mut sequence = vec![b'A'; 1600];
        sequence[1495..1504].copy_from_slice(context);
        let mut contents = b">21\n".to_vec();
        contents.extend_from_slice(&sequence);
        contents.push(b'\n');
        std::fs::write(&fa, contents).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t1600\t4\t1600\t1601\n").unwrap();
        let mut p = LofteePlugin::new();
        p.set_reference_fasta(Some(Arc::new(
            vep_fasta::IndexedFasta::from_path(&fa).unwrap(),
        )));
        p.nagnag_site(variant, "21", &coding_ctx(strand))
    }

    /// [`nagnag_probe_of`] for an acceptor SNV at 1500.
    fn nagnag_probe(context: &[u8; 9], strand: i8) -> bool {
        let alt = if context[4] == b'A' { b'G' } else { b'A' };
        let v = InputVariant::new("21".into(), 1500, 1500, vec![context[4]], vec![alt]);
        nagnag_probe_of(context, strand, &v)
    }

    /// LoF.pm:497 tests the length of the variant's span expanded by four bases
    /// each side, nine for every one-base reference span, so Ensembl VEP 116.2
    /// flags the one-base deletion of the acceptor's A (`21 25709431 TA T`, the
    /// feature `A/-` at 25709432 in the context `TTGTAGAAG`) as it flags the
    /// SNV there, while an insertion's span (`start = end + 1`) expands to eight
    /// bases and a two-base deletion's to ten, and neither is flagged. A span
    /// within four bases of the chromosome start has no nine-base context.
    #[test]
    fn nagnag_site_tests_every_one_base_span_and_no_other() {
        let context = b"TTGTAGAAG";
        let deletion = InputVariant::new("21".into(), 1500, 1500, b"A".to_vec(), b"-".to_vec());
        assert!(nagnag_probe_of(context, 1, &deletion));
        let insertion = InputVariant::new("21".into(), 1501, 1500, b"-".to_vec(), b"T".to_vec());
        assert!(!nagnag_probe_of(context, 1, &insertion));
        let two_bases = InputVariant::new("21".into(), 1500, 1501, b"AG".to_vec(), b"-".to_vec());
        assert!(!nagnag_probe_of(context, 1, &two_bases));
        let at_the_start = InputVariant::new("21".into(), 2, 2, b"A".to_vec(), b"-".to_vec());
        assert!(!nagnag_probe_of(context, 1, &at_the_start));
    }

    /// Ensembl VEP 116.2 flags an acceptor SNV whose 9-mer carries `AG.AG` at
    /// any offset (LoF.pm:498, `m/AG.AG/` unanchored): here the variant is the
    /// G of the acceptor and the rescuing `AAG` lies downstream, so the motif
    /// begins at offset 3.
    #[test]
    fn nagnag_site_matches_a_tandem_acceptor_downstream_of_the_variant() {
        assert!(nagnag_probe(b"TGTAGAAGA", 1));
    }

    /// Ensembl VEP 116.2 flags the A of the acceptor `AG` when the rescuing
    /// acceptor follows it (LoF.pm:498): the motif begins at offset 4, the last
    /// window of the 9-mer.
    #[test]
    fn nagnag_site_matches_the_acceptor_a_with_the_motif_ending_the_9mer() {
        assert!(nagnag_probe(b"TTGTAGAAG", 1));
    }

    /// Ensembl VEP 116.2 flags the A of the acceptor `AG` when the rescuing
    /// acceptor precedes it (LoF.pm:498): the motif begins at offset 1.
    #[test]
    fn nagnag_site_matches_the_acceptor_a_with_an_upstream_tandem_acceptor() {
        assert!(nagnag_probe(b"CAGCAGGTT", 1));
    }

    /// Ensembl VEP 116.2 orients the 9-mer to the transcript before the match
    /// (LoF.pm:491): the genome's `TCTTCTACA` reads `TGTAGAAGA` on a reverse
    /// strand transcript and is flagged, while the same bases on the forward
    /// strand are not.
    #[test]
    fn nagnag_site_matches_the_motif_on_the_transcript_strand() {
        assert!(nagnag_probe(b"TCTTCTACA", -1));
        assert!(!nagnag_probe(b"TCTTCTACA", 1));
    }

    /// A 9-mer with an `AG` but no second `AG` three bases on carries no
    /// `AG.AG` window, so Ensembl VEP 116.2 writes no flag (LoF.pm:498).
    #[test]
    fn nagnag_site_needs_two_acceptor_dinucleotides_three_bases_apart() {
        assert!(!nagnag_probe(b"TGTAGTTCA", 1));
    }
}
