// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Resolved configuration for a VEP run.
//!
//! Takes parsed [`Args`] and applies defaults, expands
//! shorthand flags (e.g., `--everything`), and resolves directory paths.

use crate::args::Args;
use crate::transcript_index::TranscriptIndexImpl;

/// Fully resolved configuration for a VEP run.
#[derive(Debug)]
#[allow(dead_code)] // Many fields are parsed from CLI args but only consumed by specific code paths
pub struct Config {
    pub input_file: String,
    pub input_data: Option<String>,
    pub format: String,

    pub output_file: String,
    pub output_format: String,
    /// Only consulted when `output_format == "parquet"`.
    /// `"nested"` = one row per variant, CSQ subfields as LIST<T>.
    /// `"flat"` = one row per (variant × consequence) tuple.
    pub parquet_shape: String,
    /// Keep the Parquet intermediate TSV after DuckDB finalize (debug aid).
    pub keep_intermediate: bool,
    pub parquet_row_group_size: u64,
    pub compress_output: Option<String>,
    pub force_overwrite: bool,
    pub no_headers: bool,

    pub cache: bool,
    pub offline: bool,
    pub dir: String,
    pub dir_cache: String,
    pub dir_plugins: Option<String>,
    pub cache_version: u32,
    pub fasta: Option<String>,
    pub shift_3prime: bool,

    pub species: String,
    pub assembly: Option<String>,

    pub buffer_size: usize,
    /// Structural variants wider than this many bases keep their VCF line without
    /// consequences and are dropped from the JSON output, as VEP's `--max_sv_size` does.
    pub max_sv_size: u64,
    pub fork: usize,
    /// BGZF decompression worker threads for gzipped VCF inputs (1 = single-threaded).
    pub decompression_threads: usize,
    pub distance: (u64, u64),
    /// Transcript-overlap index implementation (A/B benchmarking flag).
    /// Default [`TranscriptIndexImpl::Bin`].
    pub transcript_index_impl: TranscriptIndexImpl,

    pub symbol: bool,
    pub biotype: bool,
    pub canonical: bool,
    pub mane: bool,
    pub mane_select: bool,
    pub tsl: bool,
    pub appris: bool,
    pub ccds: bool,
    pub protein: bool,
    pub uniprot: bool,
    pub xref_refseq: bool,
    pub hgvs: bool,
    pub hgvsg: bool,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    pub numbers: bool,
    pub domains: bool,
    pub regulatory: bool,
    pub variant_class: bool,
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub pubmed: bool,
    pub gene_phenotype: bool,
    pub mirna: bool,
    pub check_existing: bool,
    pub coding_only: bool,
    pub no_intergenic: bool,

    pub pick: bool,
    pub pick_allele: bool,
    pub per_gene: bool,
    pub pick_allele_gene: bool,
    pub most_severe: bool,
    pub summary: bool,
    pub flag_pick: bool,
    pub flag_pick_allele: bool,
    pub flag_pick_allele_gene: bool,
    pub pick_order: Option<Vec<String>>,
    pub allele_number: bool,
    pub show_ref_allele: bool,
    pub minimal: bool,
    pub total_length: bool,

    pub no_stats: bool,
    pub stats_file: Option<String>,
    pub stats_text: bool,

    pub verbose: bool,
    pub quiet: bool,

    pub plugin: Vec<String>,
    pub custom: Vec<String>,

    pub refseq: bool,
    pub merged: bool,
    pub gencode_basic: bool,

    pub vcf_info_field: String,

    pub json_cache: Option<String>,

    pub fields: Option<Vec<String>>,
    pub dont_skip: bool,
    pub allow_non_variant: bool,
}

/// The user-visible notice for `--regulatory`, which vep-rs accepts for Perl VEP
/// command-line compatibility but does not implement.
///
/// No regulatory-feature cache is loaded, so none of the seven regulatory terms of
/// VEP's vocabulary (`TFBS_ablation`, `TFBS_amplification`, `TF_binding_site_variant`,
/// `regulatory_region_ablation`, `regulatory_region_amplification`,
/// `regulatory_region_variant`, `sequence_variant`) is ever emitted; the flag's only
/// effect is the extra CSQ header fields the VCF writer adds. Returns `None` when the
/// flag was neither passed nor implied by `--everything`, so a run that never asked
/// for regulatory annotation is not warned about it.
pub(crate) fn regulatory_unimplemented_notice(
    explicit: bool,
    everything: bool,
) -> Option<&'static str> {
    if explicit {
        Some(
            "--regulatory has no effect: regulatory-feature annotation is not implemented, \
             so no TFBS_* or regulatory_region_* consequence is emitted",
        )
    } else if everything {
        Some(
            "--everything implies --regulatory, which has no effect: regulatory-feature \
             annotation is not implemented, so no TFBS_* or regulatory_region_* \
             consequence is emitted",
        )
    } else {
        None
    }
}

impl Config {
    /// Build a resolved configuration from parsed CLI arguments.
    pub fn from_args(args: Args) -> Self {
        let output_format = if args.vcf {
            "vcf".to_string()
        } else if args.json {
            "json".to_string()
        } else if args.tab {
            "tab".to_string()
        } else {
            args.output_format
        };
        let parquet_shape = args.parquet_shape;
        let keep_intermediate = args.keep_intermediate;
        let parquet_row_group_size = args.parquet_row_group_size;

        let cache = args.cache || args.offline;

        let dir = expand_tilde(&args.dir);
        let dir_cache = args
            .dir_cache
            .map(|d| expand_tilde(&d))
            .unwrap_or_else(|| dir.clone());

        let cache_version = args.cache_version.unwrap_or(vep_core::VEP_VERSION);

        let distance = parse_distance(&args.distance);

        let pick_order = args.pick_order.map(|o| {
            o.split(',')
                .map(|s| s.trim().to_string())
                .collect::<Vec<_>>()
        });

        let fields = args.fields.map(|f| {
            f.split(',')
                .map(|s| s.trim().to_string())
                .collect::<Vec<_>>()
        });

        let everything = args.everything;
        // VEP's `--vcf` switches on the symbol, biotype and exon/intron-number
        // fields, which its CSQ layout names unconditionally.
        let vcf_implied = output_format == "vcf";

        let symbol = args.symbol || everything || vcf_implied;
        let biotype = args.biotype || everything || vcf_implied;
        let canonical = args.canonical || everything;
        let mane = args.mane || everything;
        let tsl = args.tsl || everything;
        let appris = args.appris || everything;
        let ccds = args.ccds || everything;
        let protein = args.protein || everything;
        let uniprot = args.uniprot || everything;
        let hgvs = args.hgvs || everything;
        let hgvsg = args.hgvsg || everything;
        let numbers = args.numbers || everything || vcf_implied;
        let domains = args.domains || everything;
        let regulatory = args.regulatory || everything;
        if let Some(notice) = regulatory_unimplemented_notice(args.regulatory, everything) {
            tracing::warn!("{}", notice);
        }
        let variant_class = args.variant_class || everything;
        let af = args.af || everything;
        let af_1kg = args.af_1kg || everything;
        let af_gnomade = args.af_gnomade || everything;
        let af_gnomadg = args.af_gnomadg || everything;
        let max_af = args.max_af || everything;
        let pubmed = args.pubmed || everything;
        let gene_phenotype = args.gene_phenotype || everything;
        let mirna = args.mirna || everything;
        let check_existing = args.check_existing || everything;
        let allele_number = args.allele_number || everything;

        let sift = args.sift.or_else(|| {
            if everything {
                Some("b".to_string())
            } else {
                None
            }
        });
        let polyphen = args.polyphen.or_else(|| {
            if everything {
                Some("b".to_string())
            } else {
                None
            }
        });

        Config {
            input_file: args.input_file,
            input_data: args.input_data,
            format: args.format,
            output_file: args.output_file,
            output_format,
            parquet_shape,
            keep_intermediate,
            parquet_row_group_size,
            compress_output: args.compress_output,
            force_overwrite: args.force_overwrite,
            no_headers: args.no_headers,
            cache,
            offline: args.offline,
            dir,
            dir_cache,
            dir_plugins: args.dir_plugins,
            cache_version,
            fasta: args.fasta,
            shift_3prime: args.shift_3prime,
            species: args.species,
            assembly: args.assembly,
            buffer_size: args.buffer_size,
            max_sv_size: args.max_sv_size,
            fork: resolve_fork(args.fork),
            decompression_threads: resolve_decompression_threads(args.decompression_threads),
            distance,
            // clap validates the flag at parse time; a failure here means the
            // validator drifted from `TranscriptIndexImpl::from_str`, and the
            // default is safer than a crash.
            transcript_index_impl: args
                .transcript_index_impl
                .parse::<TranscriptIndexImpl>()
                .unwrap_or_default(),
            symbol,
            biotype,
            canonical,
            mane,
            mane_select: args.mane_select || everything,
            tsl,
            appris,
            ccds,
            protein,
            uniprot,
            xref_refseq: args.xref_refseq,
            hgvs,
            hgvsg,
            sift,
            polyphen,
            numbers,
            domains,
            regulatory,
            variant_class,
            af,
            af_1kg,
            af_gnomade,
            af_gnomadg,
            max_af,
            pubmed,
            gene_phenotype,
            mirna,
            check_existing,
            coding_only: args.coding_only,
            no_intergenic: args.no_intergenic,
            pick: args.pick,
            pick_allele: args.pick_allele,
            per_gene: args.per_gene,
            pick_allele_gene: args.pick_allele_gene,
            most_severe: args.most_severe,
            summary: args.summary,
            flag_pick: args.flag_pick,
            flag_pick_allele: args.flag_pick_allele,
            flag_pick_allele_gene: args.flag_pick_allele_gene,
            pick_order,
            allele_number,
            show_ref_allele: args.show_ref_allele,
            minimal: args.minimal,
            total_length: args.total_length,
            no_stats: args.no_stats,
            stats_file: args.stats_file,
            stats_text: args.stats_text,
            verbose: args.verbose,
            quiet: args.quiet,
            plugin: args.plugin,
            custom: args.custom,
            refseq: args.refseq,
            merged: args.merged,
            gencode_basic: args.gencode_basic,
            vcf_info_field: args.vcf_info_field,
            json_cache: args.json_cache,
            fields,
            dont_skip: args.dont_skip,
            allow_non_variant: args.allow_non_variant,
        }
    }

    /// The cache directory the output headers name: the JSON cache when one is
    /// given, else VEP's `<dir_cache>/<species>[_<assembly>]/<version>` layout.
    pub fn cache_dir_path(&self) -> String {
        if let Some(json_cache) = &self.json_cache {
            return json_cache.clone();
        }
        let assembly_suffix = self
            .assembly
            .as_deref()
            .map(|a| format!("_{}", a))
            .unwrap_or_default();
        format!(
            "{}/{}{}/{}",
            self.dir_cache, self.species, assembly_suffix, self.cache_version
        )
    }
}

/// Auto-detect cap for `--decompression_threads=0`: decompression throughput
/// plateaus between 4 and 10 threads and is flat past that, so the auto default
/// stops at 10 rather than oversubscribing large hosts. User-supplied values are
/// not capped; any positive integer pins the worker count.
const DECOMPRESSION_THREADS_AUTO_CAP: usize = 10;

/// Resolve the --fork value. `0` means "auto": use the number of available
/// logical CPUs (no cap: annotation workers are CPU-bound and benefit from
/// the full core count on large hosts). Any non-zero value is returned
/// as-is so users can still pin the worker count explicitly.
///
/// Public so an embedding consumer of this crate can apply the same
/// auto-detect default without going through the full `Config::from_args`
/// pipeline.
pub fn resolve_fork(raw: usize) -> usize {
    if raw != 0 {
        return raw;
    }
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// Resolve the --decompression_threads value. `0` means "auto": use the
/// number of available logical CPUs, capped at
/// [`DECOMPRESSION_THREADS_AUTO_CAP`]. Any non-zero value is returned as-is
/// so users can tune the worker count for their workload and hardware.
fn resolve_decompression_threads(raw: usize) -> usize {
    if raw != 0 {
        return raw;
    }
    std::thread::available_parallelism()
        .map(|n| n.get().min(DECOMPRESSION_THREADS_AUTO_CAP))
        .unwrap_or(4)
}

/// Parse the --distance flag value into (upstream, downstream) bp.
///
/// Accepts either a single number ("5000") or a pair ("5000,2000").
fn parse_distance(s: &str) -> (u64, u64) {
    if let Some((up, down)) = s.split_once(',') {
        let up = up
            .trim()
            .parse::<u64>()
            .unwrap_or(vep_core::DEFAULT_DISTANCE);
        let down = down
            .trim()
            .parse::<u64>()
            .unwrap_or(vep_core::DEFAULT_DISTANCE);
        (up, down)
    } else {
        let d = s
            .trim()
            .parse::<u64>()
            .unwrap_or(vep_core::DEFAULT_DISTANCE);
        (d, d)
    }
}

/// Expand a leading `~` to the user's home directory.
fn expand_tilde(path: &str) -> String {
    if let Some(rest) = path.strip_prefix("~/") {
        if let Ok(home) = std::env::var("HOME") {
            return format!("{}/{}", home, rest);
        }
    }
    if path == "~" {
        if let Ok(home) = std::env::var("HOME") {
            return home;
        }
    }
    path.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;

    #[test]
    fn test_everything_expands_flags() {
        let args = Args::parse_from(["vep", "--everything"]);
        let config = Config::from_args(args);

        assert!(config.symbol);
        assert!(config.biotype);
        assert!(config.canonical);
        assert!(config.protein);
        assert!(config.hgvs);
        assert!(config.hgvsg);
        assert!(config.numbers);
        assert!(config.domains);
        assert!(config.regulatory);
        assert!(config.variant_class);
        assert!(config.af);
        assert!(config.af_1kg);
        assert!(config.af_gnomade);
        assert!(config.af_gnomadg);
        assert!(config.max_af);
        assert!(config.pubmed);
        assert!(config.uniprot);
        assert!(config.mane);
        assert!(config.mane_select);
        assert!(config.tsl);
        assert!(config.appris);
        assert!(config.gene_phenotype);
        assert!(config.mirna);
        assert!(config.ccds);
        assert!(config.check_existing);
        assert!(config.allele_number);
        assert_eq!(config.sift.as_deref(), Some("b"));
        assert_eq!(config.polyphen.as_deref(), Some("b"));
    }

    #[test]
    fn test_offline_implies_cache() {
        let args = Args::parse_from(["vep", "--offline"]);
        let config = Config::from_args(args);
        assert!(config.cache);
        assert!(config.offline);
    }

    #[test]
    fn test_vcf_shorthand() {
        let args = Args::parse_from(["vep", "--vcf"]);
        let config = Config::from_args(args);
        assert_eq!(config.output_format, "vcf");
        // VEP's --vcf implies the symbol, biotype and numbers fields.
        assert!(config.symbol && config.biotype && config.numbers);
        let plain = Config::from_args(Args::parse_from(["vep"]));
        assert!(!plain.symbol && !plain.biotype && !plain.numbers);
    }

    #[test]
    fn test_json_shorthand() {
        let args = Args::parse_from(["vep", "--json"]);
        let config = Config::from_args(args);
        assert_eq!(config.output_format, "json");
    }

    #[test]
    fn test_tab_shorthand() {
        let args = Args::parse_from(["vep", "--tab"]);
        let config = Config::from_args(args);
        assert_eq!(config.output_format, "tab");
    }

    #[test]
    fn test_parse_distance_single() {
        assert_eq!(parse_distance("5000"), (5000, 5000));
    }

    #[test]
    fn test_parse_distance_pair() {
        assert_eq!(parse_distance("5000,2000"), (5000, 2000));
    }

    #[test]
    fn test_dir_cache_defaults_to_dir() {
        let args = Args::parse_from(["vep", "--dir", "/data/vep"]);
        let config = Config::from_args(args);
        assert_eq!(config.dir_cache, "/data/vep");
    }

    #[test]
    fn test_dir_cache_explicit() {
        let args = Args::parse_from(["vep", "--dir", "/data/vep", "--dir_cache", "/cache/vep"]);
        let config = Config::from_args(args);
        assert_eq!(config.dir_cache, "/cache/vep");
    }

    #[test]
    fn test_cache_dir_path() {
        let args = Args::parse_from(["vep", "--dir_cache", "/opt/vep", "--assembly", "GRCh38"]);
        let config = Config::from_args(args);
        assert_eq!(
            config.cache_dir_path(),
            format!("/opt/vep/homo_sapiens_GRCh38/{}", vep_core::VEP_VERSION)
        );
    }

    #[test]
    fn test_expand_tilde() {
        assert_eq!(expand_tilde("/absolute/path"), "/absolute/path");
        assert_eq!(expand_tilde("relative/path"), "relative/path");
        if std::env::var("HOME").is_ok() {
            let expanded = expand_tilde("~/test");
            assert!(!expanded.starts_with('~'));
            assert!(expanded.ends_with("/test"));
        }
    }

    #[test]
    fn test_sift_explicit_overrides_everything() {
        let args = Args::parse_from(["vep", "--everything", "--sift", "p"]);
        let config = Config::from_args(args);
        assert_eq!(config.sift.as_deref(), Some("p"));
    }

    #[test]
    fn regulatory_notice_fires_only_when_requested() {
        assert!(regulatory_unimplemented_notice(false, false).is_none());
        assert!(regulatory_unimplemented_notice(true, false)
            .unwrap()
            .starts_with("--regulatory has no effect"));
        assert!(regulatory_unimplemented_notice(false, true)
            .unwrap()
            .starts_with("--everything implies --regulatory"));
    }
}
