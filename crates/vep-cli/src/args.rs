// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CLI argument parsing for VEP using clap 4 derive API.
//!
//! All flags and aliases match the Perl VEP CLI interface so that users
//! can switch between implementations with minimal changes. The Perl VEP
//! uses underscores in flag names (e.g., `--input_file`, `--dir_cache`),
//! so every underscore-containing flag has an explicit `long = "..."`.

use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "vep",
    version = env!("CARGO_PKG_VERSION"),
    about = "Ensembl Variant Effect Predictor (Rust implementation)",
    long_about = "Predict functional effects of genomic variants.\n\n\
                  Usage: vep -i input.vcf -o output.txt --json_cache <dir>"
)]
pub struct Args {
    /// Input file path (use "STDIN" for standard input)
    #[arg(short = 'i', long = "input_file", default_value = "STDIN")]
    pub input_file: String,

    /// Input data string (accepted for VEP compatibility; not implemented, input is read from --input_file)
    #[arg(long = "input_data")]
    pub input_data: Option<String>,

    /// Input file format (vcf, ensembl, hgvs, region; "guess" = auto-detect)
    #[arg(long, default_value = "guess")]
    pub format: String,

    /// Output file path (use "STDOUT" for standard output)
    #[arg(
        short = 'o',
        long = "output_file",
        default_value = "variant_effect_output.txt"
    )]
    pub output_file: String,

    /// Output format: vep, vcf, json, tab, parquet.
    /// `parquet` writes a TSV intermediate and post-processes to Parquet via a
    /// `duckdb` CLI subprocess; the `duckdb` binary must be on PATH.
    #[arg(long = "output_format", default_value = "vep")]
    pub output_format: String,

    /// Parquet row shape (only effective when --output_format=parquet).
    ///   nested (default): one row per variant, CSQ subfields are LIST<T>.
    ///   flat: one row per (variant × transcript-consequence) tuple, CSQ
    ///   subfields are scalar.
    #[arg(
        long = "parquet_shape",
        default_value = "nested",
        value_parser = ["nested", "flat"]
    )]
    pub parquet_shape: String,

    /// Keep the Parquet intermediate TSV after DuckDB finalize succeeds.
    /// Debug aid: normally the intermediate is deleted once the Parquet
    /// directory is written.
    #[arg(long = "keep_intermediate", default_value_t = false, hide = true)]
    pub keep_intermediate: bool,

    /// Rows per Parquet row group (only effective when --output_format=parquet).
    /// Rounded up to a multiple of 2,048, DuckDB's vector size, which is also
    /// the dictionary size limit. Test aid for producing several row groups
    /// from a small input.
    #[arg(
        long = "parquet_row_group_size",
        default_value_t = 122_880,
        hide = true
    )]
    pub parquet_row_group_size: u64,

    /// Shorthand for --output_format vcf
    #[arg(long)]
    pub vcf: bool,

    /// Shorthand for --output_format json
    #[arg(long)]
    pub json: bool,

    /// Shorthand for --output_format tab
    #[arg(long)]
    pub tab: bool,

    /// Compress output file (accepted for VEP compatibility; not implemented, output is written uncompressed)
    #[arg(long = "compress_output")]
    pub compress_output: Option<String>,

    /// Overwrite output file if it already exists
    #[arg(long = "force_overwrite", visible_alias = "force")]
    pub force_overwrite: bool,

    /// Suppress header lines in output
    #[arg(long = "no_headers")]
    pub no_headers: bool,

    /// Accepted for VEP compatibility; no effect (the annotation cache is --json_cache)
    #[arg(long)]
    pub cache: bool,

    /// Accepted for VEP compatibility; sets --cache, which has no effect (vep-rs always runs offline)
    #[arg(long)]
    pub offline: bool,

    /// Top-level VEP data directory (accepted for VEP compatibility; nothing is loaded from it, it only heads the cache path in the output header)
    #[arg(long, default_value = "~/.vep")]
    pub dir: String,

    /// Cache directory (defaults to --dir; accepted for VEP compatibility, nothing is loaded from it, it only names the cache path in the output header)
    #[arg(long = "dir_cache")]
    pub dir_cache: Option<String>,

    /// Plugin directory
    #[arg(long = "dir_plugins")]
    pub dir_plugins: Option<String>,

    /// Cache version named in the output header's cache path (defaults to the VEP release; selects no cache)
    #[arg(long = "cache_version")]
    pub cache_version: Option<u32>,

    /// FASTA file for reference sequence lookup
    #[arg(long, visible_alias = "fa")]
    pub fasta: Option<String>,

    /// Enable 3' shifting of insertions/deletions (requires --fasta)
    #[arg(long = "shift_3prime")]
    pub shift_3prime: bool,

    /// Species named in the output header's cache path (the annotated data is the --json_cache contents)
    #[arg(short = 's', long, default_value = "homo_sapiens")]
    pub species: String,

    /// Genome assembly (e.g., GRCh38)
    #[arg(short = 'a', long)]
    pub assembly: Option<String>,

    /// Number of variants to process per batch
    #[arg(long = "buffer_size", default_value_t = 5000)]
    pub buffer_size: usize,

    /// Structural variants spanning more than this many bases are written to the
    /// VCF output without consequences and omitted from the JSON output
    #[arg(long = "max_sv_size", default_value_t = 10_000_000)]
    pub max_sv_size: u64,

    /// Number of parallel annotation worker threads.
    /// Default 0 = auto: uses the number of available logical CPUs. Set to 1
    /// to force single-threaded annotation, or any positive value to pin the
    /// worker count.
    #[arg(long, visible_alias = "threads", default_value_t = 0)]
    pub fork: usize,

    /// Number of BGZF decompression worker threads for .gz/.bgz VCF inputs.
    /// Default 0 = auto: uses the number of available logical CPUs, capped
    /// at 10 (measured sweet-spot: throughput plateaus around 4-10 threads
    /// and is flat past that). User-supplied positive
    /// values are not capped: pass any integer to pin the worker count,
    /// e.g. `--decompression_threads=1` for single-threaded.
    #[arg(long = "decompression_threads", default_value_t = 0)]
    pub decompression_threads: usize,

    /// Transcript overlap index implementation (hidden, for A/B benchmarking).
    /// Options: bin (UCSC binning, default), coitrees (implicit interval tree),
    /// sorted (sorted array + suffix_max_end). All three must produce identical
    /// overlap sets; any difference between implementations indicates a bug.
    #[arg(
        long = "transcript_index_impl",
        default_value = "bin",
        hide = true,
        value_parser = ["bin", "coitrees", "sorted"]
    )]
    pub transcript_index_impl: String,

    /// Upstream/downstream distance in bp (or "up,down" for asymmetric)
    #[arg(long, default_value = "5000")]
    pub distance: String,

    /// Shorthand to enable commonly used output fields
    #[arg(short = 'e', long)]
    pub everything: bool,

    /// Add gene symbol to output
    #[arg(long)]
    pub symbol: bool,

    /// Add transcript biotype to output
    #[arg(long)]
    pub biotype: bool,

    /// Flag canonical transcripts
    #[arg(long)]
    pub canonical: bool,

    /// Add MANE Select and MANE Plus Clinical flags
    #[arg(long)]
    pub mane: bool,

    /// Add MANE Select flag only
    #[arg(long = "mane_select")]
    pub mane_select: bool,

    /// Add Transcript Support Level
    #[arg(long)]
    pub tsl: bool,

    /// Add APPRIS isoform annotation
    #[arg(long)]
    pub appris: bool,

    /// Add CCDS transcript ID
    #[arg(long)]
    pub ccds: bool,

    /// Add Ensembl protein ID
    #[arg(long)]
    pub protein: bool,

    /// Add UniProt cross-references
    #[arg(long)]
    pub uniprot: bool,

    /// Add RefSeq transcript cross-references
    #[arg(long = "xref_refseq")]
    pub xref_refseq: bool,

    /// Add HGVS coding and protein notations
    #[arg(long)]
    pub hgvs: bool,

    /// Add the HGVSg column (accepted for VEP compatibility; the value is not computed and stays empty)
    #[arg(long)]
    pub hgvsg: bool,

    /// SIFT predictions (p = prediction, s = score, b = both)
    #[arg(long)]
    pub sift: Option<String>,

    /// PolyPhen predictions (p = prediction, s = score, b = both)
    #[arg(long)]
    pub polyphen: Option<String>,

    /// Add exon/intron numbers
    #[arg(long)]
    pub numbers: bool,

    /// Add overlapping protein domain information
    #[arg(long)]
    pub domains: bool,

    /// Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns
    #[arg(long)]
    pub regulatory: bool,

    /// Add variant class from SO
    #[arg(long = "variant_class")]
    pub variant_class: bool,

    /// Add global allele frequency from 1000 Genomes Phase 3
    #[arg(long)]
    pub af: bool,

    /// Add continental allele frequencies from 1000 Genomes
    #[arg(long = "af_1kg")]
    pub af_1kg: bool,

    /// Add gnomAD exome allele frequencies
    #[arg(long = "af_gnomade")]
    pub af_gnomade: bool,

    /// Add gnomAD genome allele frequencies
    #[arg(long = "af_gnomadg")]
    pub af_gnomadg: bool,

    /// Add highest allele frequency across populations
    #[arg(long = "max_af")]
    pub max_af: bool,

    /// Add PubMed IDs for co-located variants
    #[arg(long)]
    pub pubmed: bool,

    /// Add gene-phenotype associations
    #[arg(long = "gene_phenotype")]
    pub gene_phenotype: bool,

    /// Add the miRNA column (accepted for VEP compatibility; the value is not computed and stays empty)
    #[arg(long)]
    pub mirna: bool,

    /// Check for co-located known variants
    #[arg(long = "check_existing")]
    pub check_existing: bool,

    /// Accepted for VEP compatibility; not implemented
    #[arg(long = "coding_only")]
    pub coding_only: bool,

    /// Accepted for VEP compatibility; not implemented
    #[arg(long = "no_intergenic")]
    pub no_intergenic: bool,

    /// Pick one consequence per variant
    #[arg(long)]
    pub pick: bool,

    /// Pick one consequence per variant allele
    #[arg(long = "pick_allele")]
    pub pick_allele: bool,

    /// Pick one consequence per gene
    #[arg(long = "per_gene")]
    pub per_gene: bool,

    /// Pick one consequence per allele per gene
    #[arg(long = "pick_allele_gene")]
    pub pick_allele_gene: bool,

    /// Report only the single most severe consequence
    #[arg(long = "most_severe")]
    pub most_severe: bool,

    /// Report a summary of all consequences
    #[arg(long)]
    pub summary: bool,

    /// Flag (but don't remove) non-picked consequences
    #[arg(long = "flag_pick")]
    pub flag_pick: bool,

    /// Flag non-picked consequences per allele
    #[arg(long = "flag_pick_allele")]
    pub flag_pick_allele: bool,

    /// Flag non-picked consequences per allele per gene
    #[arg(long = "flag_pick_allele_gene")]
    pub flag_pick_allele_gene: bool,

    /// Accepted for VEP compatibility; not implemented (the pick order is fixed to VEP's default)
    #[arg(long = "pick_order")]
    pub pick_order: Option<String>,

    /// Add allele number from input VCF
    #[arg(long = "allele_number")]
    pub allele_number: bool,

    /// Show reference allele in output
    #[arg(long = "show_ref_allele")]
    pub show_ref_allele: bool,

    /// Convert alleles to minimal representation
    #[arg(long)]
    pub minimal: bool,

    /// Accepted for VEP compatibility; not implemented
    #[arg(long = "total_length")]
    pub total_length: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it)
    #[arg(long = "no_stats")]
    pub no_stats: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(long = "stats_file")]
    pub stats_file: Option<String>,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(long = "stats_text")]
    pub stats_text: bool,

    /// Verbose output (debug-level logging)
    #[arg(short = 'v', long)]
    pub verbose: bool,

    /// Suppress status/progress messages
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Load plugin(s) - can be specified multiple times
    #[arg(long)]
    pub plugin: Vec<String>,

    /// Custom annotation source(s) (accepted for VEP compatibility; not implemented, sources are not loaded)
    #[arg(long)]
    pub custom: Vec<String>,

    /// Accepted for VEP compatibility; not implemented (the transcript set is the --json_cache contents)
    #[arg(long)]
    pub refseq: bool,

    /// Accepted for VEP compatibility; not implemented (the transcript set is the --json_cache contents)
    #[arg(long)]
    pub merged: bool,

    /// Accepted for VEP compatibility; not implemented (the transcript set is the --json_cache contents)
    #[arg(long = "gencode_basic")]
    pub gencode_basic: bool,

    /// VCF INFO field name for VEP annotation
    #[arg(long = "vcf_info_field", default_value = "CSQ")]
    pub vcf_info_field: String,

    /// Path to the JSON transcript cache directory (built by vep-cache-builder or converted from a Perl VEP cache)
    #[arg(long = "json_cache")]
    pub json_cache: Option<String>,

    /// Accepted for VEP compatibility; not implemented (the columns follow the field flags)
    #[arg(long)]
    pub fields: Option<String>,

    /// Abort the run on a variant that fails parsing instead of skipping it (VEP keeps the variant)
    #[arg(long = "dont_skip")]
    pub dont_skip: bool,

    /// Process non-variant lines in VCF (e.g., reference-only calls)
    #[arg(long = "allow_non_variant")]
    pub allow_non_variant: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_args() {
        let args = Args::parse_from(["vep"]);
        assert_eq!(args.input_file, "STDIN");
        assert_eq!(args.output_file, "variant_effect_output.txt");
        assert_eq!(args.output_format, "vep");
        assert_eq!(args.format, "guess");
        assert_eq!(args.species, "homo_sapiens");
        assert_eq!(args.buffer_size, 5000);
        assert_eq!(args.fork, 0);
        assert!(!args.cache);
        assert!(!args.offline);
        assert!(!args.everything);
    }

    #[test]
    fn test_input_output_flags() {
        let args = Args::parse_from([
            "vep",
            "-i",
            "test.vcf",
            "-o",
            "out.txt",
            "--cache",
            "--offline",
        ]);
        assert_eq!(args.input_file, "test.vcf");
        assert_eq!(args.output_file, "out.txt");
        assert!(args.cache);
        assert!(args.offline);
    }

    #[test]
    fn test_force_alias() {
        let args = Args::parse_from(["vep", "--force"]);
        assert!(args.force_overwrite);
    }

    #[test]
    fn test_threads_alias() {
        let args = Args::parse_from(["vep", "--threads", "4"]);
        assert_eq!(args.fork, 4);
    }

    #[test]
    fn test_fa_alias() {
        let args = Args::parse_from(["vep", "--fa", "/path/to/ref.fa"]);
        assert_eq!(args.fasta.as_deref(), Some("/path/to/ref.fa"));
    }

    #[test]
    fn test_everything_flag() {
        let args = Args::parse_from(["vep", "--everything"]);
        assert!(args.everything);
    }

    #[test]
    fn test_output_format_shortcuts() {
        let args = Args::parse_from(["vep", "--vcf"]);
        assert!(args.vcf);

        let args = Args::parse_from(["vep", "--json"]);
        assert!(args.json);

        let args = Args::parse_from(["vep", "--tab"]);
        assert!(args.tab);
    }

    #[test]
    fn test_repeated_plugin_flags() {
        let args = Args::parse_from([
            "vep",
            "--plugin",
            "CADD,/path/to/cadd.tsv.gz",
            "--plugin",
            "dbNSFP,/path/to/dbnsfp.gz",
        ]);
        assert_eq!(args.plugin.len(), 2);
    }

    #[test]
    fn test_underscore_flag_names() {
        let args = Args::parse_from([
            "vep",
            "--dir_cache",
            "/cache",
            "--dir_plugins",
            "/plugins",
            "--cache_version",
            "115",
            "--buffer_size",
            "1000",
            "--output_format",
            "json",
            "--force_overwrite",
            "--no_headers",
        ]);
        assert_eq!(args.dir_cache.as_deref(), Some("/cache"));
        assert_eq!(args.dir_plugins.as_deref(), Some("/plugins"));
        assert_eq!(args.cache_version, Some(115));
        assert_eq!(args.buffer_size, 1000);
        assert_eq!(args.output_format, "json");
        assert!(args.force_overwrite);
        assert!(args.no_headers);
    }

    #[test]
    fn test_output_control_underscore_flags() {
        let args = Args::parse_from([
            "vep",
            "--mane_select",
            "--xref_refseq",
            "--variant_class",
            "--af_1kg",
            "--af_gnomade",
            "--af_gnomadg",
            "--max_af",
            "--gene_phenotype",
            "--check_existing",
            "--coding_only",
            "--no_intergenic",
        ]);
        assert!(args.mane_select);
        assert!(args.xref_refseq);
        assert!(args.variant_class);
        assert!(args.af_1kg);
        assert!(args.af_gnomade);
        assert!(args.af_gnomadg);
        assert!(args.max_af);
        assert!(args.gene_phenotype);
        assert!(args.check_existing);
        assert!(args.coding_only);
        assert!(args.no_intergenic);
    }

    #[test]
    fn test_filtering_underscore_flags() {
        let args = Args::parse_from([
            "vep",
            "--pick_allele",
            "--per_gene",
            "--pick_allele_gene",
            "--most_severe",
            "--flag_pick",
            "--flag_pick_allele",
            "--flag_pick_allele_gene",
            "--pick_order",
            "canonical,biotype,rank",
            "--allele_number",
            "--show_ref_allele",
            "--total_length",
            "--dont_skip",
            "--allow_non_variant",
        ]);
        assert!(args.pick_allele);
        assert!(args.per_gene);
        assert!(args.pick_allele_gene);
        assert!(args.most_severe);
        assert!(args.flag_pick);
        assert!(args.flag_pick_allele);
        assert!(args.flag_pick_allele_gene);
        assert_eq!(args.pick_order.as_deref(), Some("canonical,biotype,rank"));
        assert!(args.allele_number);
        assert!(args.show_ref_allele);
        assert!(args.total_length);
        assert!(args.dont_skip);
        assert!(args.allow_non_variant);
    }
}
