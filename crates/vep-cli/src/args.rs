// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CLI argument parsing for VEP using clap 4 derive API.
//!
//! All flags and aliases match the Perl VEP CLI interface so that users
//! can switch between implementations with minimal changes. The Perl VEP
//! uses underscores in flag names (e.g., `--input_file`, `--dir_cache`),
//! so every underscore-containing flag has an explicit `long = "..."`.
//!
//! A flag the error table of [`crate::config`] governs (refused, or accepted
//! without effect) carries no help text here: [`Args::command_with_error_table`]
//! applies the table's, so `--help` and the generated CLI reference read one
//! source. Every flag lives under a help heading, the section of the reference
//! it is documented in.
//!
//! The command line reaches this parser through [`crate::argv::prepass`], which
//! adds the flags of the `--config` file, `vep.ini` and the `VEP_*` variables
//! and spells single-dash and other-case long options as declared here; a
//! unique prefix of a long name is accepted here as well (`infer_long_args`),
//! as Getopt::Long's `auto_abbrev` accepts it.

use clap::{ArgAction, CommandFactory, Parser};

/// Ensembl VEP's `--max_sv_size` default: the span above which a structural variant
/// is oversize, and the repeat-sequence length above which a tandem repeat stays
/// symbolic.
pub const DEFAULT_MAX_SV_SIZE: i64 = vep_core::sv_record::DEFAULT_MAX_SV_SIZE as i64;

/// The help's first line is the binary's name and version, as `--version` prints.
const HELP_TEMPLATE: &str =
    "{name} {version}\n{about-with-newline}\n{usage-heading} {usage}\n\n{all-args}{after-help}";

#[derive(Parser, Debug)]
#[command(
    name = "vep",
    version = env!("CARGO_PKG_VERSION"),
    about = "Ensembl Variant Effect Predictor (Rust implementation)",
    long_about = "Predict functional effects of genomic variants.\n\n\
                  Usage: vep -i input.vcf -o output.txt --json_cache <dir>",
    help_template = HELP_TEMPLATE,
    args_override_self = true,
    infer_long_args = true
)]
pub struct Args {
    /// Configuration file, one flag per line as `<name> <value>` without the dashes (`format vcf`, `offline 1`); a flag the command line also gives keeps the command line's value
    #[arg(help_heading = "Input", long, value_name = "FILE")]
    pub config: Option<String>,

    /// Input file path (use "STDIN" for standard input)
    #[arg(
        help_heading = "Input",
        short = 'i',
        long = "input_file",
        default_value = "STDIN"
    )]
    pub input_file: String,

    /// Input data string (accepted for VEP compatibility; not implemented, input is read from --input_file)
    #[arg(help_heading = "Input", long = "input_data")]
    pub input_data: Option<String>,

    /// Input file format ("guess" = auto-detect)
    #[arg(
        help_heading = "Input",
        long,
        default_value = "guess",
        value_parser = ["ensembl", "vcf", "hgvs", "id", "spdi", "region", "guess"]
    )]
    pub format: String,

    #[arg(help_heading = "Input", long)]
    pub gp: bool,

    #[arg(help_heading = "Input", long = "ambiguous_hgvs")]
    pub ambiguous_hgvs: bool,

    #[arg(help_heading = "Input", long)]
    pub phased: bool,

    #[arg(help_heading = "Input", long = "no_check_variants_order")]
    pub no_check_variants_order: bool,

    /// Output file path (use "STDOUT" for standard output, in any case); the last of several wins
    #[arg(
        help_heading = "Output",
        short = 'o',
        long = "output_file",
        action = ArgAction::Append,
        default_value = "variant_effect_output.txt"
    )]
    pub output_file: Vec<String>,

    /// Output format, matched in any case; `ensembl` is another name for the default `vep` format. `parquet` is vep-rs's own: it writes a TSV intermediate and post-processes it
    /// to Parquet through a `duckdb` CLI subprocess; the `duckdb` binary must be on PATH
    #[arg(
        help_heading = "Output",
        long = "output_format",
        default_value = "vep",
        ignore_case = true,
        value_parser = ["vep", "ensembl", "vcf", "json", "tab", "parquet"]
    )]
    pub output_format: String,

    /// Parquet row shape (only effective when --output_format=parquet).
    ///   nested (default): one row per variant, CSQ subfields are LIST<T>.
    ///   flat: one row per (variant × transcript-consequence) tuple, CSQ
    ///   subfields are scalar.
    #[arg(
        help_heading = "Output",
        long = "parquet_shape",
        default_value = "nested",
        value_parser = ["nested", "flat"]
    )]
    pub parquet_shape: String,

    /// Keep the Parquet intermediate TSV after DuckDB finalize succeeds.
    /// Debug aid: normally the intermediate is deleted once the Parquet
    /// directory is written.
    #[arg(
        help_heading = "Output",
        long = "keep_intermediate",
        default_value_t = false,
        hide = true
    )]
    pub keep_intermediate: bool,

    /// Rows per Parquet row group (only effective when --output_format=parquet).
    /// Rounded up to a multiple of 2,048, DuckDB's vector size, which is also
    /// the dictionary size limit. Test aid for producing several row groups
    /// from a small input.
    #[arg(
        help_heading = "Output",
        long = "parquet_row_group_size",
        default_value_t = 122_880,
        hide = true
    )]
    pub parquet_row_group_size: u64,

    /// Shorthand for --output_format vcf; also switches on --symbol, --biotype and --numbers
    #[arg(help_heading = "Output", long)]
    pub vcf: bool,

    /// Shorthand for --output_format json
    #[arg(help_heading = "Output", long)]
    pub json: bool,

    /// Shorthand for --output_format tab
    #[arg(help_heading = "Output", long)]
    pub tab: bool,

    /// Compress the output, STDOUT included: `gzip` writes a single-member gzip stream, `bgzip` a BGZF stream (a gzip file that bgzip, tabix and samtools read); the output file keeps the name given. Not accepted with --output_format parquet, which carries its own compression
    #[arg(help_heading = "Output", long = "compress_output", value_parser = ["gzip", "bgzip"])]
    pub compress_output: Option<String>,

    /// Overwrite output file if it already exists
    #[arg(
        help_heading = "Output",
        long = "force_overwrite",
        visible_alias = "force"
    )]
    pub force_overwrite: bool,

    /// Suppress header lines in output
    #[arg(help_heading = "Output", long = "no_headers")]
    pub no_headers: bool,

    /// Print the cache directory in the output headers as [PATH]/<leaf>
    #[arg(help_heading = "Output", long = "mask_header_cache_path")]
    pub mask_header_cache_path: bool,

    /// File the warnings of the run are written to, created at the first warning;
    /// default <output_file>_warnings.txt beside the output (STDOUT_warnings.txt in the
    /// working directory for -o STDOUT). STDERR, in any case, keeps them on standard
    /// error alone. A warning is also echoed on standard error unless --quiet
    #[arg(help_heading = "Output", long = "warning_file", value_name = "FILE")]
    pub warning_file: Option<String>,

    /// File listing every input record the run skipped: a header line naming the
    /// input, then one line per record with its line number, the record's first 40
    /// characters and the reason. An existing file is refused before the run without
    /// --force_overwrite; the file is written when the run ends
    #[arg(
        help_heading = "Output",
        long = "skipped_variants_file",
        value_name = "FILE"
    )]
    pub skipped_variants_file: Option<String>,

    /// The columns of the --tab output and the subfields of the VCF CSQ entry, comma-separated, in this order and replacing the whole default list; any field name a flag or a plugin produces is valid, and a name nothing produces prints `-` (tab) or an empty subfield (VCF) under a `?` description. The default and JSON formats ignore it; --output_format parquet refuses it
    #[arg(help_heading = "Output", long)]
    pub fields: Option<String>,

    /// VCF INFO field name for VEP annotation
    #[arg(
        help_heading = "Output",
        long = "vcf_info_field",
        default_value = "CSQ"
    )]
    pub vcf_info_field: String,

    /// Keep an input VCF's existing annotation: its `##INFO` and `##VEP` header lines and the INFO entry under --vcf_info_field's name stay, and the new entry is appended under the same name, so the line carries both
    #[arg(help_heading = "Output", long = "keep_csq")]
    pub keep_csq: bool,

    #[arg(help_heading = "Output", long = "keep_ann")]
    pub keep_ann: bool,

    #[arg(help_heading = "Output", long)]
    pub solr: bool,

    #[arg(help_heading = "Output", long, value_name = "FORMAT")]
    pub convert: Option<String>,

    /// Annotate from the local cache, which vep-rs always does; a directory value is the cache root, as --dir_cache, when that flag is absent
    #[arg(
        help_heading = "Cache / Offline",
        long,
        num_args = 0..=1,
        default_missing_value = crate::config::CACHE_FLAG_ALONE,
        value_name = "DIR"
    )]
    pub cache: Option<String>,

    #[arg(help_heading = "Cache / Offline", long)]
    pub offline: bool,

    /// Top-level VEP data directory: the cache root when neither --dir_cache nor a --cache directory is given
    #[arg(help_heading = "Cache / Offline", long, default_value = "~/.vep")]
    pub dir: String,

    /// Cache root: the cache read is <dir_cache>/<species>[_refseq|_merged]/<cache_version>_<assembly>, as Ensembl VEP lays it out (defaults to --dir)
    #[arg(help_heading = "Cache / Offline", long = "dir_cache")]
    pub dir_cache: Option<String>,

    /// Cache version: the prefix of the cache directory under the species directory (defaults to the VEP release)
    #[arg(help_heading = "Cache / Offline", long = "cache_version")]
    pub cache_version: Option<u32>,

    /// The cache directory itself, the one holding info.json and transcripts/ (built by vep-cache-builder or converted from a Perl VEP cache); it is read in place of the --dir_cache layout, and only --assembly is checked against its info.json
    #[arg(
        help_heading = "Cache / Offline",
        long = "json_cache",
        visible_alias = "full_cache_dir"
    )]
    pub json_cache: Option<String>,

    /// Print the cache's source versions as sorted key<TAB>value lines and exit
    #[arg(help_heading = "Cache / Offline", long = "show_cache_info")]
    pub show_cache_info: bool,

    /// FASTA file for reference sequence lookup; without it, the first .fa by name in the cache directory that has its .fai beside it is used
    #[arg(help_heading = "Cache / Offline", long, visible_alias = "fa")]
    pub fasta: Option<String>,

    /// Do not use a FASTA file found in the cache directory
    #[arg(help_heading = "Cache / Offline", long = "no_fasta")]
    pub no_fasta: bool,

    /// Move each insertion or deletion to its most 3' position on the transcript, over
    /// a 1,000-base window of the reference, before its consequence terms, positions,
    /// codons and amino acids are computed (1), or annotate it where the input placed it
    /// (0); needs a FASTA and stops the run without one. The HGVS notations shift under
    /// --shift_hgvs, not this flag. The shift runs along the genome only: a transcript
    /// whose cached sequence is edited away from the genome is not shifted again along
    /// its own sequence
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_3prime",
        value_parser = ["0", "1"]
    )]
    pub shift_3prime: Option<String>,

    /// Move each insertion or deletion to its most 3' position on the transcript before
    /// its HGVSc and HGVSp are built (1, the default), or describe it where the input
    /// placed it (0), which also sets --shift_3prime to 0 and cannot be combined with
    /// --shift_3prime or --shift_genomic
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_hgvs",
        value_parser = ["0", "1"]
    )]
    pub shift_hgvs: Option<String>,

    /// Under --shift_3prime 1, print the shifted position in Location and the shifted
    /// sequence in Allele on every row (1), intergenic rows moving towards higher
    /// coordinates, or keep the input position (0); without --shift_3prime 1 nothing
    /// moves
    #[arg(
        help_heading = "Cache / Offline",
        long = "shift_genomic",
        value_parser = ["0", "1"]
    )]
    pub shift_genomic: Option<String>,

    /// Add the SHIFT_LENGTH column: under --shift_3prime 1, the bases each insertion
    /// or deletion moved along the transcript, 0 on every other transcript row
    #[arg(help_heading = "Cache / Offline", long = "shift_length")]
    pub shift_length: bool,

    #[arg(help_heading = "Cache / Offline", long)]
    pub sereal: bool,

    /// Species: names the cache directory under --dir_cache; the Latin name only, with no alias such as human
    #[arg(
        help_heading = "Species / Assembly",
        short = 's',
        long,
        default_value = "homo_sapiens"
    )]
    pub species: String,

    /// Genome assembly (e.g., GRCh38): chooses among several caches of one version and must equal the cache's own name exactly; hg19 and hg38 are read as GRCh37 and GRCh38
    #[arg(help_heading = "Species / Assembly", short = 'a', long)]
    pub assembly: Option<String>,

    /// Annotate against GRCh37: sets --assembly GRCh37, over an explicit --assembly, as Ensembl VEP's option set does
    #[arg(help_heading = "Species / Assembly", long)]
    pub grch37: bool,

    #[arg(help_heading = "Database", long)]
    pub database: bool,

    #[arg(help_heading = "Database", long, value_name = "HOST")]
    pub host: Option<String>,

    #[arg(help_heading = "Database", long, value_name = "PORT")]
    pub port: Option<String>,

    #[arg(help_heading = "Database", short = 'u', long, value_name = "USER")]
    pub user: Option<String>,

    #[arg(
        help_heading = "Database",
        long,
        visible_alias = "pass",
        value_name = "PASSWORD"
    )]
    pub password: Option<String>,

    #[arg(help_heading = "Database", long, value_name = "FILE")]
    pub registry: Option<String>,

    #[arg(help_heading = "Database", long = "db_version", value_name = "VERSION")]
    pub db_version: Option<String>,

    #[arg(help_heading = "Database", long)]
    pub genomes: bool,

    #[arg(
        help_heading = "Database",
        long = "is_multispecies",
        value_name = "0|1"
    )]
    pub is_multispecies: Option<String>,

    /// Number of variants to process per batch
    #[arg(help_heading = "Runtime", long = "buffer_size", default_value_t = 5000)]
    pub buffer_size: usize,

    /// A structural variant spanning more than this many bases is oversize; -1 lifts
    /// the limit, as in Ensembl VEP 116. A <CNV:TR> record whose repeat sequence is
    /// longer than this stays symbolic rather than reading a direction of change
    /// (under -1 every one does), as in Ensembl VEP 116. Ensembl VEP 116 drops an
    /// oversize record before annotation and writes nothing for it in any format;
    /// vep-rs annotates it: rows in the default, tab and parquet outputs, its VCF
    /// line without consequences, no JSON object (docs/output-formats.md and
    /// docs/intended-divergences.md describe the divergence)
    #[arg(
        help_heading = "Runtime",
        long = "max_sv_size",
        default_value_t = DEFAULT_MAX_SV_SIZE,
        allow_negative_numbers = true,
        value_parser = clap::value_parser!(i64).range(-1..)
    )]
    pub max_sv_size: i64,

    /// Number of parallel annotation worker threads (Ensembl VEP forks processes; vep-rs runs threads and is parallel by default).
    /// Default 0 = auto: uses the number of available logical CPUs, capped at
    /// 32 (wall time stops improving past that and rises on wider hosts). Set
    /// to 1 to force single-threaded annotation, or any positive value to pin
    /// the worker count; user-supplied values are not capped.
    #[arg(
        help_heading = "Runtime",
        long,
        visible_alias = "threads",
        default_value_t = 0
    )]
    pub fork: usize,

    /// Number of BGZF decompression worker threads for .gz/.bgz VCF inputs.
    /// Default 0 = auto: uses the number of available logical CPUs, capped
    /// at 10 (measured sweet-spot: throughput plateaus around 4-10 threads
    /// and is flat past that). User-supplied positive
    /// values are not capped: pass any integer to pin the worker count,
    /// e.g. `--decompression_threads=1` for single-threaded.
    #[arg(
        help_heading = "Runtime",
        long = "decompression_threads",
        default_value_t = 0
    )]
    pub decompression_threads: usize,

    /// Transcript overlap index implementation (hidden, for A/B benchmarking).
    /// Options: bin (UCSC binning, default), coitrees (implicit interval tree),
    /// sorted (sorted array + suffix_max_end). All three must produce identical
    /// overlap sets; any difference between implementations indicates a bug.
    #[arg(
        help_heading = "Runtime",
        long = "transcript_index_impl",
        default_value = "bin",
        hide = true,
        value_parser = ["bin", "coitrees", "sorted"]
    )]
    pub transcript_index_impl: String,

    /// Upstream and downstream distance in bp: one integer, or the pair "up,down"
    /// (a downstream part written 0 reads as the upstream value)
    #[arg(help_heading = "Runtime", long, default_value = "5000")]
    pub distance: String,

    /// Shorthand for Ensembl VEP's --everything set: --sift b, --polyphen b, --ccds, --hgvs, --symbol,
    /// --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg,
    /// --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris,
    /// --variant_class, --gene_phenotype and --mirna, each overriding an explicit value as in Ensembl VEP.
    /// Every column of Ensembl VEP's --everything except the regulatory columns: --regulatory is
    /// switched on and annotates nothing, so no RegulatoryFeature or MotifFeature row is written
    #[arg(help_heading = "Output Control", short = 'e', long)]
    pub everything: bool,

    /// Add gene symbol to output
    #[arg(help_heading = "Output Control", long)]
    pub symbol: bool,

    /// Add transcript biotype to output
    #[arg(help_heading = "Output Control", long)]
    pub biotype: bool,

    /// Flag canonical transcripts
    #[arg(help_heading = "Output Control", long)]
    pub canonical: bool,

    /// Add MANE Select and MANE Plus Clinical flags
    #[arg(help_heading = "Output Control", long)]
    pub mane: bool,

    /// Add MANE Select flag only
    #[arg(help_heading = "Output Control", long = "mane_select")]
    pub mane_select: bool,

    /// Add Transcript Support Level
    #[arg(help_heading = "Output Control", long)]
    pub tsl: bool,

    /// Add APPRIS isoform annotation
    #[arg(help_heading = "Output Control", long)]
    pub appris: bool,

    /// Add CCDS transcript ID
    #[arg(help_heading = "Output Control", long)]
    pub ccds: bool,

    /// Add Ensembl protein ID
    #[arg(help_heading = "Output Control", long)]
    pub protein: bool,

    /// Add the translation version to the Ensembl protein ID (ENSP00000286808.3); needs --protein
    #[arg(help_heading = "Output Control", long = "protein_version")]
    pub protein_version: bool,

    /// Add UniProt cross-references
    #[arg(help_heading = "Output Control", long)]
    pub uniprot: bool,

    /// Add RefSeq transcript cross-references
    #[arg(help_heading = "Output Control", long = "xref_refseq")]
    pub xref_refseq: bool,

    /// Add HGVS coding and protein notations; needs a FASTA named with --fasta or found in the cache directory and stops the run without one
    #[arg(help_heading = "Output Control", long)]
    pub hgvs: bool,

    /// Add the HGVS genomic notation (HGVSg); needs a FASTA named with --fasta or found in the cache directory and stops the run without one
    #[arg(help_heading = "Output Control", long)]
    pub hgvsg: bool,

    /// SIFT predictions (p = prediction, s = score, b = both)
    #[arg(help_heading = "Output Control", long, value_parser = ["s", "p", "b"])]
    pub sift: Option<String>,

    /// PolyPhen predictions (p = prediction, s = score, b = both)
    #[arg(help_heading = "Output Control", long, value_parser = ["s", "p", "b"])]
    pub polyphen: Option<String>,

    #[arg(help_heading = "Output Control", long, value_name = "CONDEL")]
    pub condel: Option<String>,

    /// Add exon/intron numbers
    #[arg(help_heading = "Output Control", long)]
    pub numbers: bool,

    /// Add overlapping protein domain information
    #[arg(help_heading = "Output Control", long)]
    pub domains: bool,

    /// Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns
    #[arg(help_heading = "Output Control", long)]
    pub regulatory: bool,

    #[arg(
        help_heading = "Output Control",
        long = "regulatory_gff",
        value_name = "FILE"
    )]
    pub regulatory_gff: Option<String>,

    #[arg(help_heading = "Output Control", long = "extended_promoters")]
    pub extended_promoters: bool,

    /// Add variant class from SO
    #[arg(help_heading = "Output Control", long = "variant_class")]
    pub variant_class: bool,

    /// Add the OverlapBP and OverlapPC columns to the tab, VCF and Parquet outputs: the
    /// bases of a structural variant inside the transcript and the percentage of the
    /// transcript they cover (transcript rows only; regulatory features are not annotated)
    #[arg(help_heading = "Output Control", long)]
    pub overlaps: bool,

    /// Check for co-located known variants
    #[arg(help_heading = "Output Control", long = "check_existing")]
    pub check_existing: bool,

    /// Report a known variant as co-located when its start and end equal the input variant's, whatever its alleles; such a variant prints on every allele's row
    #[arg(help_heading = "Output Control", long = "no_check_alleles")]
    pub no_check_alleles: bool,

    /// Never report a known variant whose alleles the cache does not hold (COSMIC and HGMD records, dbSNP records without variation); without the flag such a variant is reported when its start and end equal the input variant's
    #[arg(help_heading = "Output Control", long = "exclude_null_alleles")]
    pub exclude_null_alleles: bool,

    /// Report known variants whose failed value is at most this: 0 excludes variants that failed Ensembl's quality checks, 1 includes them
    #[arg(
        help_heading = "Output Control",
        long,
        value_name = "INT",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    pub failed: i64,

    /// Accepted for VEP compatibility; not implemented
    #[arg(help_heading = "Output Control", long = "coding_only")]
    pub coding_only: bool,

    /// Accepted for VEP compatibility; not implemented
    #[arg(help_heading = "Output Control", long = "no_intergenic")]
    pub no_intergenic: bool,

    #[arg(help_heading = "Output Control", long)]
    pub hgnc: bool,

    /// Add global allele frequency from 1000 Genomes Phase 3
    #[arg(help_heading = "Allele Frequencies", long)]
    pub af: bool,

    /// Add continental allele frequencies from 1000 Genomes
    #[arg(help_heading = "Allele Frequencies", long = "af_1kg")]
    pub af_1kg: bool,

    /// Add gnomAD exome allele frequencies
    #[arg(help_heading = "Allele Frequencies", long = "af_gnomade")]
    pub af_gnomade: bool,

    /// Add gnomAD genome allele frequencies
    #[arg(help_heading = "Allele Frequencies", long = "af_gnomadg")]
    pub af_gnomadg: bool,

    /// Add highest allele frequency across populations
    #[arg(help_heading = "Allele Frequencies", long = "max_af")]
    pub max_af: bool,

    /// Add PubMed IDs for co-located variants
    #[arg(help_heading = "Allele Frequencies", long)]
    pub pubmed: bool,

    /// Rewrite the AFR, AMR, ASN and EUR frequency values of each co-located variant as a bare minor-allele frequency: the allele is stripped, the first value kept and a value above 0.5 inverted; the allele-keyed frequency columns then find no allele in them and print empty
    #[arg(help_heading = "Allele Frequencies", long = "old_maf")]
    pub old_maf: bool,

    /// Add gene-phenotype associations
    #[arg(help_heading = "Allele Frequencies", long = "gene_phenotype")]
    pub gene_phenotype: bool,

    /// Pick one consequence per variant
    #[arg(help_heading = "Filtering / Picking", long)]
    pub pick: bool,

    /// Pick one consequence per variant allele
    #[arg(help_heading = "Filtering / Picking", long = "pick_allele")]
    pub pick_allele: bool,

    /// Pick one consequence per gene
    #[arg(help_heading = "Filtering / Picking", long = "per_gene")]
    pub per_gene: bool,

    /// Pick one consequence per allele per gene
    #[arg(help_heading = "Filtering / Picking", long = "pick_allele_gene")]
    pub pick_allele_gene: bool,

    /// Report only the single most severe consequence
    #[arg(help_heading = "Filtering / Picking", long = "most_severe")]
    pub most_severe: bool,

    /// Report a summary of all consequences
    #[arg(help_heading = "Filtering / Picking", long)]
    pub summary: bool,

    /// Flag (but don't remove) non-picked consequences
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick")]
    pub flag_pick: bool,

    /// Flag non-picked consequences per allele
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick_allele")]
    pub flag_pick_allele: bool,

    /// Flag non-picked consequences per allele per gene
    #[arg(help_heading = "Filtering / Picking", long = "flag_pick_allele_gene")]
    pub flag_pick_allele_gene: bool,

    /// The criteria the pick family compares, in order, as a comma-separated list; a
    /// shorter list compares its criteria and then keeps the first transcript in
    /// stable-id order. Default: mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq
    #[arg(
        help_heading = "Filtering / Picking",
        long = "pick_order",
        value_delimiter = ',',
        value_parser = [
            "mane_select",
            "mane_plus_clinical",
            "canonical",
            "appris",
            "tsl",
            "biotype",
            "ccds",
            "rank",
            "length",
            "ensembl",
            "refseq",
        ]
    )]
    pub pick_order: Option<Vec<String>>,

    /// Annotate only the transcripts that pass a filter, in Ensembl VEP's filter
    /// grammar: units of `field operator value` joined by `and` and `or`, grouped
    /// by parentheses and inverted by `not`; repeated flags are joined by `and`.
    /// The fields are the cache's transcript keys, named exactly, in another
    /// case, without the leading `_`, or by a unique prefix: `stable_id`,
    /// `version`, `biotype`, `start`, `end`, `strand`, `is_canonical`, `source`,
    /// `dbID`, `cdna_coding_start`, `cdna_coding_end`, `_gene_stable_id`,
    /// `_gene_symbol`, `_gene_symbol_source`, `_gene_hgnc_id`, `_ccds`,
    /// `_refseq`, `_protein`, `_swissprot`, `_trembl`, `_uniparc`,
    /// `_gene_phenotype`; the object keys `translation` (coding transcripts),
    /// `slice`, `attributes`, `_gene`, `_trans_exon_array`,
    /// `_variation_effect_feature_cache` and `_vep_lazy_loaded`, which are defined
    /// and read meaningfully by `exists` and `nex` alone; and `description`,
    /// `coding_region_start`, `coding_region_end`, `_gene_version` and
    /// `_uniprot_isoform`, which are undefined (the cache leaves the first three
    /// undefined, the transcript loader does not read the gene version, and the
    /// converted cache does not carry the isoform list).
    /// `_source_cache` is absent, as it is from every cache but a merged one. The
    /// operators are `is` (`eq`, `=`), `ne` (`!=`), `>`, `<`, `>=`, `<=`, `exists`
    /// (`defined`, the default), `nex`, `match` (`matches`, `regex`, `re`), `nre`,
    /// `in` (a comma-separated list or a file with one value per line) and
    /// `is_child`. `match` and `nre` take the `regex` crate's dialect,
    /// case-insensitive and unanchored: a back-reference or an escape the crate
    /// lacks is an error naming the pattern. A record whose every transcript is
    /// filtered out is intergenic.
    #[arg(
        help_heading = "Filtering / Picking",
        long = "transcript_filter",
        value_name = "FILTER"
    )]
    pub transcript_filter: Vec<String>,

    /// Add allele number from input VCF
    #[arg(help_heading = "Filtering / Picking", long = "allele_number")]
    pub allele_number: bool,

    /// Show reference allele in output
    #[arg(help_heading = "Filtering / Picking", long = "show_ref_allele")]
    pub show_ref_allele: bool,

    /// Convert alleles to minimal representation; switches on --allele_number as in Ensembl VEP
    #[arg(help_heading = "Filtering / Picking", long)]
    pub minimal: bool,

    /// Accepted for VEP compatibility; not implemented
    #[arg(help_heading = "Filtering / Picking", long = "total_length")]
    pub total_length: bool,

    /// Load plugin(s) - can be specified multiple times
    #[arg(help_heading = "Plugins", long)]
    pub plugin: Vec<String>,

    /// Custom annotation source(s) (accepted for VEP compatibility; not implemented, sources are not loaded)
    #[arg(help_heading = "Plugins", long)]
    pub custom: Vec<String>,

    #[arg(help_heading = "Plugins", long = "custom_suppress_filter")]
    pub custom_suppress_filter: bool,

    /// Directory the plugin libraries are looked up in; default <--dir>/Plugins. A
    /// directory that does not exist is warned about when --plugin names a plugin
    /// (an error under --safe)
    #[arg(help_heading = "Plugins", long = "dir_plugins", value_name = "DIR")]
    pub dir_plugins: Option<String>,

    /// Make a plugin that cannot be loaded, or a missing plugins directory, an error
    /// instead of a warning, and skip the check of incompatible flag pairs
    #[arg(help_heading = "Plugins", long)]
    pub safe: bool,

    /// Reads the cache under the species directory suffixed _refseq (the transcripts annotated are that cache's contents); no effect under --json_cache
    #[arg(help_heading = "Transcript Sets", long)]
    pub refseq: bool,

    /// Reads the cache under the species directory suffixed _merged (the transcripts annotated are that cache's contents); no effect under --json_cache
    #[arg(help_heading = "Transcript Sets", long)]
    pub merged: bool,

    /// Annotate against the GENCODE basic transcript set alone: a transcript without a
    /// gencode_basic attribute is dropped as the cache loads, before the transcript index
    /// is built
    #[arg(help_heading = "Transcript Sets", long = "gencode_basic")]
    pub gencode_basic: bool,

    /// Annotate against the GENCODE primary transcript set alone: a transcript without a
    /// gencode_primary attribute is dropped as the cache loads; the set exists for human
    /// GRCh38 only, so --assembly naming another assembly is warned about
    #[arg(help_heading = "Transcript Sets", long = "gencode_primary")]
    pub gencode_primary: bool,

    /// Add the GENCODE_PRIMARY column: 1 on a transcript of the GENCODE primary set, empty
    /// otherwise; the set exists for human GRCh38 only, so --assembly naming another
    /// assembly is warned about
    #[arg(help_heading = "Transcript Sets", long = "flag_gencode_primary")]
    pub flag_gencode_primary: bool,

    #[arg(help_heading = "Transcript Sets", long = "use_given_ref")]
    pub use_given_ref: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it)
    #[arg(help_heading = "Statistics", long = "no_stats")]
    pub no_stats: bool,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(help_heading = "Statistics", long = "stats_file", visible_alias = "sf")]
    pub stats_file: Option<String>,

    /// Accepted for VEP compatibility; not implemented (no statistics file is written)
    #[arg(help_heading = "Statistics", long = "stats_text")]
    pub stats_text: bool,

    /// Verbose output (debug-level logging)
    #[arg(help_heading = "Verbosity", short = 'v', long)]
    pub verbose: bool,

    /// Suppress the status and progress lines and the standard-error copy of each
    /// warning; the warnings file is still written
    #[arg(help_heading = "Verbosity", short = 'q', long)]
    pub quiet: bool,

    /// Abort the run on a variant that fails parsing instead of skipping it (VEP keeps the variant)
    #[arg(help_heading = "Miscellaneous", long = "dont_skip")]
    pub dont_skip: bool,

    /// Process non-variant lines in VCF (e.g., reference-only calls)
    #[arg(help_heading = "Miscellaneous", long = "allow_non_variant")]
    pub allow_non_variant: bool,

    /// Add the miRNA column (accepted for VEP compatibility; the value is not computed and stays empty)
    #[arg(help_heading = "Miscellaneous", long)]
    pub mirna: bool,

    #[arg(help_heading = "Miscellaneous", long)]
    pub lrg: bool,

    #[arg(help_heading = "Miscellaneous", long, value_name = "DIR")]
    pub tmpdir: Option<String>,
}

impl Args {
    /// The command with the error table's help text applied to every flag the
    /// table governs; the parser `main` runs and the source of `--help` and the
    /// CLI reference.
    pub fn command_with_error_table() -> clap::Command {
        Args::command().mut_args(|arg| {
            let entry = arg
                .get_long()
                .and_then(crate::config::error_table_entry)
                .filter(|e| e.declared);
            match entry {
                Some(entry) => arg.help(entry.help()),
                None => arg,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{Config, ParseFailure};

    #[test]
    fn test_default_args() {
        let args = Args::parse_from(["vep"]);
        assert_eq!(args.input_file, "STDIN");
        assert_eq!(args.output_file, ["variant_effect_output.txt"]);
        assert_eq!(args.output_format, "vep");
        assert_eq!(args.format, "guess");
        assert_eq!(args.species, "homo_sapiens");
        assert_eq!(args.buffer_size, 5000);
        assert_eq!(args.fork, 0);
        assert!(args.cache.is_none());
        assert!(!args.offline);
        assert!(!args.everything);
    }

    /// `--cache` takes an optional value (`cache:s`, Config.pm:238): alone it is
    /// on, and a following path is its directory.
    #[test]
    fn cache_takes_an_optional_directory() {
        let args = Args::parse_from(["vep", "--cache", "--offline"]);
        assert_eq!(args.cache.as_deref(), Some("1"));
        assert!(args.offline);
        let args = Args::parse_from(["vep", "--cache", "/opt/vep", "--offline"]);
        assert_eq!(args.cache.as_deref(), Some("/opt/vep"));
        let args = Args::parse_from(["vep", "--cache=/opt/vep"]);
        assert_eq!(args.cache.as_deref(), Some("/opt/vep"));
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
        assert_eq!(args.output_file, ["out.txt"]);
        assert_eq!(args.cache.as_deref(), Some("1"));
        assert!(args.offline);
    }

    /// A repeated single-value flag keeps its last value, as Getopt::Long does;
    /// `-o` accumulates (`output_file|o=s@`) so the duplicate can be reported.
    #[test]
    fn repeated_flags_keep_the_last_value() {
        let args = Args::parse_from(["vep", "-i", "a.vcf", "-i", "b.vcf", "-o", "x", "-o", "y"]);
        assert_eq!(args.input_file, "b.vcf");
        assert_eq!(args.output_file, ["x", "y"]);
        let args = Args::parse_from([
            "vep", "--fork", "2", "--fork", "3", "--sift", "p", "--sift", "b",
        ]);
        assert_eq!(args.fork, 3);
        assert_eq!(args.sift.as_deref(), Some("b"));
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

    /// `-u` and `--pass` are Ensembl VEP's short forms of `--user` and `--password`
    /// (`Config.pm:90-91`), so they reach the same refusal.
    #[test]
    fn database_short_forms_parse_as_their_flags() {
        let args = Args::parse_from(["vep", "-u", "anonymous", "--pass", "secret"]);
        assert_eq!(args.user.as_deref(), Some("anonymous"));
        assert_eq!(args.password.as_deref(), Some("secret"));
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

    /// `--plugin` accumulates (`plugin=s@`, Config.pm:271) and the plugins load in
    /// command-line order (Runner.pm:925 iterates the list as given).
    #[test]
    fn repeated_plugin_flags_load_in_order() {
        let args = Args::parse_from([
            "vep",
            "--plugin",
            "CADD,/path/to/cadd.tsv.gz",
            "--plugin",
            "dbNSFP,/path/to/dbnsfp.gz",
        ]);
        assert_eq!(
            args.plugin,
            ["CADD,/path/to/cadd.tsv.gz", "dbNSFP,/path/to/dbnsfp.gz"]
        );
    }

    /// The warning sink's flags parse: `--safe` (Config.pm:272), `--warning_file`
    /// (`:173`), `--skipped_variants_file` (`:174`).
    #[test]
    fn warning_sink_flags_parse() {
        let args = Args::parse_from([
            "vep",
            "--safe",
            "--warning_file",
            "w.txt",
            "--skipped_variants_file",
            "s.txt",
        ]);
        assert!(args.safe);
        assert_eq!(args.warning_file.as_deref(), Some("w.txt"));
        assert_eq!(args.skipped_variants_file.as_deref(), Some("s.txt"));
        let args = Args::parse_from(["vep"]);
        assert!(!args.safe);
        assert!(args.warning_file.is_none());
        assert!(args.skipped_variants_file.is_none());
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

    /// `--failed` is Ensembl VEP's `failed=i` (Config.pm:119, default 0 at :308): an
    /// integer, negative included, and a syntax error for anything else; the two
    /// co-location switches parse beside it.
    #[test]
    fn failed_takes_an_integer() {
        let config = Config::from_argv(["vep", "--check_existing"]).unwrap();
        assert_eq!(config.failed, 0);
        assert!(!config.no_check_alleles);
        assert!(!config.exclude_null_alleles);
        let config = Config::from_argv([
            "vep",
            "--check_existing",
            "--failed",
            "1",
            "--no_check_alleles",
            "--exclude_null_alleles",
        ])
        .unwrap();
        assert_eq!(config.failed, 1);
        assert!(config.no_check_alleles);
        assert!(config.exclude_null_alleles);
        assert_eq!(
            Config::from_argv(["vep", "--failed=-1"]).unwrap().failed,
            -1
        );
        assert_eq!(
            Config::from_argv(["vep", "--failed", "-1"]).unwrap().failed,
            -1
        );
        let err = Config::from_argv(["vep", "--failed", "yes"]).unwrap_err();
        assert!(
            matches!(err, ParseFailure::Usage(_)),
            "a non-integer is clap's own error: {err}"
        );
    }

    #[test]
    fn protein_version_and_mask_header_cache_path_parse() {
        let config =
            Config::from_argv(["vep", "--protein_version", "--mask_header_cache_path"]).unwrap();
        assert!(config.protein_version);
        assert!(config.mask_header_cache_path);
    }

    /// `--regulatory_gff` is refused: regulatory annotation is not supported.
    #[test]
    fn regulatory_gff_is_refused() {
        let err = Config::from_argv(["vep", "--regulatory_gff", "reg.gff3.gz"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --regulatory_gff is not supported: it would annotate regulatory features \
             and motifs from a GFF3 file in place of the cache, and regulatory annotation \
             is not supported"
        );
    }

    /// `--extended_promoters` is refused: regulatory annotation is not supported.
    #[test]
    fn extended_promoters_is_refused() {
        let err = Config::from_argv(["vep", "--extended_promoters"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --extended_promoters is not supported: it would widen the promoters read \
             with --regulatory_gff to their extended bounds, and regulatory annotation is \
             not supported"
        );
    }

    /// `--custom_suppress_filter` is refused: custom annotation is not supported.
    #[test]
    fn custom_suppress_filter_is_refused() {
        let err = Config::from_argv(["vep", "--custom_suppress_filter"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [cli] --custom_suppress_filter is not supported: it would keep a --custom VCF \
             source's FILTER column out of the output unless the --custom source's fields \
             argument names it, and custom annotation is not supported"
        );
    }

    /// `--transcript_filter` repeats, each value a whole filter string with its
    /// spaces, as Ensembl VEP's `transcript_filter=s@` takes it.
    #[test]
    fn transcript_filter_repeats() {
        let args = Args::parse_from([
            "vep",
            "--transcript_filter",
            "biotype is protein_coding",
            "--transcript_filter",
            "not stable_id match ENST0000040",
        ]);
        assert_eq!(
            args.transcript_filter,
            [
                "biotype is protein_coding",
                "not stable_id match ENST0000040"
            ]
        );
        assert!(Args::parse_from(["vep"]).transcript_filter.is_empty());
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
        assert_eq!(
            args.pick_order.as_deref(),
            Some(
                &[
                    "canonical".to_string(),
                    "biotype".to_string(),
                    "rank".to_string()
                ][..]
            )
        );
        assert!(args.allele_number);
        assert!(args.show_ref_allele);
        assert!(args.total_length);
        assert!(args.dont_skip);
        assert!(args.allow_non_variant);
    }

    /// `shift_hgvs=i` (Config.pm:201) takes a value, and only 0 or 1 is read as a
    /// setting (BaseRunner.pm:498): the bare flag and any other value are syntax
    /// errors, so neither can pass as a silent default.
    #[test]
    fn shift_hgvs_takes_0_or_1() {
        assert_eq!(
            Args::parse_from(["vep", "--shift_hgvs", "0"])
                .shift_hgvs
                .as_deref(),
            Some("0")
        );
        assert_eq!(
            Args::parse_from(["vep", "--shift_hgvs=1"])
                .shift_hgvs
                .as_deref(),
            Some("1")
        );
        assert_eq!(Args::parse_from(["vep"]).shift_hgvs, None);
        for argv in [
            &["vep", "--shift_hgvs"][..],
            &["vep", "--shift_hgvs", "2"][..],
            &["vep", "--shift_hgvs", "yes"][..],
        ] {
            let err = Args::try_parse_from(argv).unwrap_err();
            assert_eq!(err.exit_code(), 2, "{argv:?}: {err}");
        }
    }

    /// `shift_3prime=i` and `shift_genomic=i` (Config.pm:231-232) take a value as
    /// `--shift_hgvs` does; the bare `--shift_3prime` is a syntax error, and
    /// `--shift_length` is a plain switch (Config.pm:233).
    #[test]
    fn shift_3prime_takes_0_or_1() {
        let args = Args::parse_from(["vep", "--shift_3prime", "1", "--shift_genomic", "0"]);
        assert_eq!(args.shift_3prime.as_deref(), Some("1"));
        assert_eq!(args.shift_genomic.as_deref(), Some("0"));
        assert!(!args.shift_length);
        assert_eq!(
            Args::parse_from(["vep", "--shift_3prime=0"])
                .shift_3prime
                .as_deref(),
            Some("0")
        );
        assert!(Args::parse_from(["vep", "--shift_length"]).shift_length);
        for argv in [
            &["vep", "--shift_3prime"][..],
            &["vep", "--shift_3prime", "2"][..],
            &["vep", "--shift_genomic"][..],
            &["vep", "--shift_genomic", "-1"][..],
        ] {
            let err = Args::try_parse_from(argv).unwrap_err();
            assert_eq!(err.exit_code(), 2, "{argv:?}: {err}");
        }
    }

    /// Every flag has help text: its own doc comment, or the error table's.
    #[test]
    fn every_flag_has_help() {
        let mut command = Args::command_with_error_table();
        command.build();
        let without: Vec<&str> = command
            .get_arguments()
            .filter(|a| a.get_help().is_none() && a.get_long_help().is_none())
            .filter_map(|a| a.get_long())
            .collect();
        assert!(without.is_empty(), "flags without help: {without:?}");
    }

    /// Every flag is documented under a help heading, one of the reference's sections.
    #[test]
    fn every_flag_has_a_help_heading() {
        let mut command = Args::command();
        command.build();
        let without: Vec<&str> = command
            .get_arguments()
            .filter(|a| a.get_long().is_some_and(|l| l != "help" && l != "version"))
            .filter(|a| a.get_help_heading().is_none())
            .filter_map(|a| a.get_long())
            .collect();
        assert!(
            without.is_empty(),
            "flags without a help heading: {without:?}"
        );
    }
}
