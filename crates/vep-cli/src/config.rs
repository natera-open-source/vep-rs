// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Resolved configuration for a VEP run.
//!
//! Takes parsed [`Args`] and applies defaults, Ensembl VEP's option sets and
//! option rules ([`crate::option_rules`]), the error table of flags vep-rs
//! refuses or accepts without effect, and resolves directory paths. Every
//! refusal is a [`CliError`], the one shape the binary prints before exiting 1.

use std::ffi::OsString;
use std::fmt;

use clap::FromArgMatches;

use crate::args::Args;
use crate::option_rules::{FlagState, OWN_RULES};
use crate::transcript_index::TranscriptIndexImpl;

/// A refusal under the error contract: printed as one line on standard error,
/// `ERROR: [<component>] <message>`, and the run exits 1. `<message>` is Ensembl
/// VEP's text where it has one, without its own `ERROR: ` prefix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CliError {
    pub component: &'static str,
    pub message: String,
}

impl CliError {
    pub fn new(component: &'static str, message: impl Into<String>) -> Self {
        CliError {
            component,
            message: message.into(),
        }
    }

    /// A command-line refusal: an unsupported flag or an option-rule violation.
    pub fn cli(message: impl Into<String>) -> Self {
        CliError::new("cli", message)
    }

    /// The line an error of the annotation run prints: a [`CliError`] as itself,
    /// any other failure under the `run` component with its cause chain.
    pub fn from_run(error: &anyhow::Error) -> String {
        match error.downcast_ref::<CliError>() {
            Some(e) => e.to_string(),
            None => CliError::new("run", format!("{error:#}")).to_string(),
        }
    }
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ERROR: [{}] {}", self.component, self.message)
    }
}

impl std::error::Error for CliError {}

/// Why a command line was not accepted: a syntax error clap reports in its own
/// words (exit 2; `--help` and `--version` arrive here too, exit 0), or a
/// refusal under the error contract (exit 1).
#[derive(Debug)]
pub enum ParseFailure {
    Usage(Box<clap::Error>),
    Refused(CliError),
}

impl fmt::Display for ParseFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseFailure::Usage(e) => write!(f, "{e}"),
            ParseFailure::Refused(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for ParseFailure {}

/// How vep-rs treats a flag the error table governs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Treatment {
    /// The flag is refused with this message under this component.
    Refused {
        component: &'static str,
        message: &'static str,
    },
    /// The flag parses and changes nothing, for this reason.
    AcceptedWithoutEffect { reason: &'static str },
}

/// One row of the error table: a flag of Ensembl VEP's command line that vep-rs
/// either refuses or accepts without effect, with what Ensembl VEP does with it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ErrorTableEntry {
    /// The long flag, without its dashes.
    pub flag: &'static str,
    /// Whether the parser declares the flag; a flag absent from Ensembl VEP
    /// 116.2's `@VEP_PARAMS` is caught by a scan of the command line before parsing.
    pub declared: bool,
    pub treatment: Treatment,
    /// What Ensembl VEP does with the flag, cited to its source at release/116.2.
    pub ensembl_vep: &'static str,
    /// Where the behaviour lives in vep-rs, or what to use instead.
    pub route: &'static str,
}

impl ErrorTableEntry {
    /// The help text `--help` and the CLI reference print for the flag.
    pub fn help(&self) -> String {
        match self.treatment {
            Treatment::Refused { message, .. } => format!("Refused with an error: {message}"),
            Treatment::AcceptedWithoutEffect { reason } => {
                format!("Accepted without effect: {reason}")
            }
        }
    }
}

const fn refused(
    flag: &'static str,
    message: &'static str,
    ensembl_vep: &'static str,
    route: &'static str,
) -> ErrorTableEntry {
    ErrorTableEntry {
        flag,
        declared: true,
        treatment: Treatment::Refused {
            component: "cli",
            message,
        },
        ensembl_vep,
        route,
    }
}

const fn noop(
    flag: &'static str,
    reason: &'static str,
    ensembl_vep: &'static str,
) -> ErrorTableEntry {
    ErrorTableEntry {
        flag,
        declared: true,
        treatment: Treatment::AcceptedWithoutEffect { reason },
        ensembl_vep,
        route: "nothing changes",
    }
}

/// The error table: every flag vep-rs refuses, declared or caught before parsing,
/// and the flags accepted without effect whose rows it owns. Help text, the
/// refusal and the CLI reference all read it.
pub const ERROR_TABLE: &[ErrorTableEntry] = &[
    refused(
        "database",
        "--database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:237 declares it; AnnotationSourceAdaptor.pm:146-200 opens the database",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "host",
        "--host is not supported: it names the Ensembl database host, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:88, default ensembldb.ensembl.org at :287",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "port",
        "--port is not supported: it names the Ensembl database port, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:89, default 3306 at :289, 3337 under --grch37 at :401",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "user",
        "--user is not supported: it names the Ensembl database user, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:90 (user|u=s), default anonymous at :288",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "password",
        "--password is not supported: it is the Ensembl database password, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:91 (password|pass=s)",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "registry",
        "--registry is not supported: it names an Ensembl registry file of database connections, and vep-rs reads a local cache only; name the cache with --json_cache, or --cache with --dir_cache",
        "Config.pm:87",
        "--json_cache, or --cache with --dir_cache",
    ),
    refused(
        "db_version",
        "--db_version is not supported: it selects an Ensembl database version; --cache_version selects a cache",
        "Config.pm:92; CacheDir.pm:323 falls back to it for the cache version",
        "--cache_version",
    ),
    refused(
        "genomes",
        "--genomes is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only",
        "Config.pm:95; the option set at :390-396 sets the host and port",
        "a local cache",
    ),
    refused(
        "is_multispecies",
        "--is_multispecies is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only",
        "Config.pm:101; BaseVEP.pm:282 reads it for the registry species",
        "a local cache",
    ),
    refused(
        "lrg",
        "--lrg is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode)",
        "Config.pm:221, incompatible with --offline at :631; Runner.pm:766 refuses it offline",
        "none offline",
    ),
    refused(
        "convert",
        "--convert has been deprecated",
        "Config.pm:636-638 lists it as deprecated and :977-985 refuses it with this text",
        "none",
    ),
    refused(
        "gp",
        "--gp is not supported: it reads 1000 Genomes pilot coordinates from INFO/GP, which no current VCF carries",
        "Config.pm:120; Parser/VCF.pm:271-285 takes chr:pos from INFO/GP and skips lines without it",
        "the record's own coordinates",
    ),
    refused(
        "ambiguous_hgvs",
        "--ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g.) only",
        "Config.pm:202; Parser/HGVS.pm:178-186, 224 resolves the protein notation through the database",
        "genomic HGVS input",
    ),
    ErrorTableEntry {
        flag: "af_exac",
        declared: false,
        treatment: Treatment::Refused {
            component: "cli",
            message: "--af_exac was removed in Ensembl VEP release 115: use --af_gnomade",
        },
        ensembl_vep: "Config.pm:618 names it among --database's incompatibilities while @VEP_PARAMS (:73-276) does not declare it, so Ensembl VEP refuses it as an unexpected parameter (vep:30)",
        route: "--af_gnomade",
    },
    refused(
        "regulatory_gff",
        "--regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported",
        "Config.pm:210; the option set at :434-439 switches on --regulatory",
        "none",
    ),
    refused(
        "extended_promoters",
        "--extended_promoters is not supported: it would widen the promoters read with --regulatory_gff to their extended bounds, and regulatory annotation is not supported",
        "Config.pm:211",
        "none",
    ),
    refused(
        "custom_suppress_filter",
        "--custom_suppress_filter is not supported: it would keep a --custom VCF source's FILTER column out of the output unless the --custom source's fields argument names it, and custom annotation is not supported",
        "Config.pm:268",
        "none",
    ),
    noop(
        "condel",
        "Ensembl VEP declares it and reads nothing from it",
        "Config.pm:206 (condel=s), read by no module",
    ),
    noop(
        "hgnc",
        "Ensembl VEP declares it and reads nothing from it",
        "Config.pm:189, read by no module",
    ),
    noop(
        "keep_ann",
        "Ensembl VEP declares it as a synonym of --keep_csq and reads nothing from it",
        "Config.pm:220, read by no module",
    ),
    noop(
        "phased",
        "Ensembl VEP marks the genotypes as phased and nothing reads the mark",
        "Config.pm:146; Parser/VCF.pm:107, 654, 734 set phased on the variant and nothing reads it",
    ),
    noop(
        "sereal",
        "the cache's serialiser is read from its own metadata",
        "Config.pm:249; CacheDir.pm:147 takes the serialiser from the cache's info file",
    ),
    noop(
        "solr",
        "Ensembl VEP declares it and writes no Solr output",
        "Config.pm:215, read by no module",
    ),
    noop(
        "tmpdir",
        "no module of Ensembl VEP reads it, and vep-rs writes its Parquet intermediate beside the output file",
        "Config.pm:259 (tmpdir=s), default /tmp at :304, read by no module",
    ),
    noop(
        "offline",
        "vep-rs always runs offline and never contacts a database",
        "Config.pm:245; the option set at :464-469 switches on --cache and off --database",
    ),
    noop(
        "use_given_ref",
        "vep-rs always uses the reference allele the input gives",
        "Config.pm:255; the option set at :588-592 switches off use_transcript_ref, which CacheDir.pm:412-416 enables for a BAM-edited cache",
    ),
    noop(
        "no_check_variants_order",
        "vep-rs never checks the input order: every record is annotated where it stands, and an input whose chromosomes interleave keeps more transcript indexes in memory",
        "Config.pm:83; InputBuffer.pm:189-258 aborts after 100 records more than 100 bp behind their predecessor",
    ),
];

/// The value `--cache` carries when given without a directory
/// ([`Args::cache`]'s `default_missing_value`).
pub const CACHE_FLAG_ALONE: &str = "1";

/// Environment variables vep-rs reads that Ensembl VEP has no counterpart for.
pub const ENVIRONMENT: &[(&str, &str)] = &[(
    "VEP_FASTA_STORAGE",
    "how the FASTA is held: `mmap` (default) maps the file and inflates a bgzipped one block by block as queries need them, `memory` reads it into RAM, a bgzipped one inflated whole; any other value is an error",
)];

/// The components an `ERROR: [<component>]` line names, with their meaning.
pub const COMPONENTS: &[(&str, &str)] = &[
    (
        "cli",
        "an unsupported flag, or a combination the option rules refuse",
    ),
    (
        "output",
        "the output file or format, the warnings file or the skipped-variants file",
    ),
    (
        "plugin",
        "under --safe, a plugin that cannot be loaded or a plugins directory that does not exist",
    ),
    (
        "cache",
        "the cache: nothing resolves at the configured path, its layout or metadata is unusable, its assembly disagrees with --assembly or with its own directory name, or a variation file cannot be opened or read",
    ),
    (
        "fasta",
        "the reference FASTA: the file or directory named cannot be read as an indexed FASTA, or a flag that reads it (`--hgvs`, `--hgvsg`, `--shift_3prime 1`) was given and none was named or found in the cache directory",
    ),
    (
        "filter",
        "a `--transcript_filter` Ensembl VEP would refuse, or whose regular expression is outside the supported dialect",
    ),
    (
        "config",
        "the `--config` file, or a `vep.ini` that is a file, cannot be read; or a configuration value the parser accepts but the run cannot use, such as a `--distance` part that is not an integer",
    ),
    (
        "run",
        "any other failure of the annotation run, with its cause chain",
    ),
];

/// The error table's entry for a flag, when it has one.
pub fn error_table_entry(flag: &str) -> Option<&'static ErrorTableEntry> {
    ERROR_TABLE.iter().find(|e| e.flag == flag)
}

/// Where a command line's warnings go: `--warning_file`, else
/// `<output_file>_warnings.txt` beside the last `-o` value as the user spelled it
/// (BaseVEP.pm:890); readable before the configuration is resolved, so the sink
/// can take the warnings the resolution itself raises.
pub fn warning_file_for(args: &Args) -> String {
    match &args.warning_file {
        Some(file) => file.clone(),
        None => crate::warnings::default_warning_file(
            args.output_file
                .last()
                .map(String::as_str)
                .unwrap_or("variant_effect_output.txt"),
        ),
    }
}

/// A parsed command line: the argument struct and the flag state clap saw.
pub struct ParsedArgs {
    pub args: Args,
    pub state: FlagState,
}

/// Parses a command line (the binary name first). A flag the error table refuses
/// before parsing is reported first, unless the line also asks for the help or
/// the version, which win and drop the flag; then clap's own syntax errors; the
/// result carries the arguments and the flag state the option rules act on.
pub fn parse_command_line<I, T>(argv: I) -> Result<ParsedArgs, ParseFailure>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let mut argv = argv.into_iter().map(Into::into);
    let program = argv.next();
    let mut options: Vec<OsString> = argv.collect();
    // The usage wins over an unexpected parameter (vep:29 before vep:30): an
    // undeclared flag of the error table is dropped from a command line that asks
    // for the help or the version, and refused otherwise. A flag outside the
    // table is clap's, which reads the line in order and refuses the flag unless
    // the help or the version comes first.
    let asks_help_or_version = options.iter().any(|a| {
        matches!(
            a.to_string_lossy().as_ref(),
            "--help" | "-h" | "--version" | "-V"
        )
    });
    for entry in ERROR_TABLE.iter().filter(|e| !e.declared) {
        let Treatment::Refused { component, message } = entry.treatment else {
            continue;
        };
        let long = format!("--{}", entry.flag);
        let with_value = format!("{long}=");
        let names_it = |a: &OsString| {
            let a = a.to_string_lossy();
            a == long || a.starts_with(&with_value)
        };
        if asks_help_or_version {
            options.retain(|a| !names_it(a));
        } else if options.iter().any(names_it) {
            return Err(ParseFailure::Refused(CliError::new(component, message)));
        }
    }
    let usage = |error: clap::Error| ParseFailure::Usage(Box::new(error));
    let mut command = Args::command_with_error_table();
    command.build();
    let matches = command
        .clone()
        .try_get_matches_from(program.into_iter().chain(options))
        .map_err(usage)?;
    let args = Args::from_arg_matches(&matches).map_err(usage)?;
    let state = FlagState::from_matches(&command, &matches);
    Ok(ParsedArgs { args, state })
}

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
    /// Print the cache directory in the output headers as `[PATH]/<leaf>`.
    pub mask_header_cache_path: bool,
    /// Where the run's warnings go: the `--warning_file` value, else
    /// `<output_file>_warnings.txt` with the output file as the user spelled it;
    /// `stderr` in any case names standard error (BaseVEP.pm:878-897).
    pub warning_file: String,
    pub skipped_variants_file: Option<String>,

    pub cache: bool,
    pub offline: bool,
    pub dir: String,
    /// The cache root the Ensembl VEP layout is resolved under: `--dir_cache`,
    /// else a directory given to `--cache`, else `--dir`.
    pub dir_cache: String,
    /// The plugins directory: `--dir_plugins`, else `<dir>/Plugins`
    /// (Runner.pm:907), a leading `~` expanded.
    pub dir_plugins: String,
    /// `--safe`: a plugin that cannot be loaded or a missing plugins directory is
    /// an error rather than a warning (Runner.pm:912-916, 929-953).
    pub safe: bool,
    pub cache_version: u32,
    /// Print the cache's source versions and exit.
    pub show_cache_info: bool,
    pub fasta: Option<String>,
    /// A directory of FASTA files to pick the species' and assembly's from.
    pub fasta_dir: Option<String>,
    /// Do not pick up a FASTA file found in the cache directory.
    pub no_fasta: bool,
    /// The `--synonyms` file, a leading `~` expanded; it replaces the cache's own
    /// `chr_synonyms.txt` (`CacheDir.pm:254` reads that file only when the flag is
    /// absent).
    pub synonyms: Option<String>,
    /// `--shift_3prime 1`: insertions and deletions are moved to their most 3'
    /// position on each transcript before annotation; off under `--shift_hgvs 0`,
    /// which sets it to 0 as `post_setup_checks` does (Runner.pm:792-795).
    pub shift_3prime: bool,
    /// `--shift_hgvs`: the HGVS notations describe the 3'-shifted allele (the
    /// default); `--shift_hgvs 0` describes the allele where the input placed it.
    pub shift_hgvs: bool,
    /// `--shift_genomic 1`: Location and Allele print the shifted position and
    /// allele under `--shift_3prime 1`.
    pub shift_genomic: bool,
    /// `--shift_length`: the SHIFT_LENGTH column is added.
    pub shift_length: bool,

    pub species: String,
    pub assembly: Option<String>,

    pub buffer_size: usize,
    /// The `--max_sv_size` limit; `None` is the flag's `-1`, no limit. Ensembl VEP
    /// 116 drops a wider structural variant before annotation and writes nothing
    /// for it in any format (`Parser.pm` `validate_vf`); vep-rs annotates it in the
    /// default, tab and Parquet outputs, carries its VCF line without consequences
    /// and omits it from the JSON output (`docs/intended-divergences.md`).
    pub max_sv_size: Option<u64>,
    pub fork: usize,
    /// BGZF decompression worker threads for gzipped VCF inputs (1 = single-threaded).
    pub decompression_threads: usize,
    pub distance: (u64, u64),
    /// Transcript-overlap index implementation (A/B benchmarking flag).
    /// Default [`TranscriptIndexImpl::Bin`].
    pub transcript_index_impl: TranscriptIndexImpl,

    /// `--everything` was given: the HGVS flags it switched on are dropped, not
    /// refused, when no FASTA is found (Runner.pm:745-754).
    pub everything: bool,
    pub symbol: bool,
    pub biotype: bool,
    pub canonical: bool,
    pub mane: bool,
    pub mane_select: bool,
    pub tsl: bool,
    pub appris: bool,
    pub ccds: bool,
    pub protein: bool,
    /// Append the translation version to the ENSP identifier (`--protein_version`).
    pub protein_version: bool,
    pub uniprot: bool,
    pub xref_refseq: bool,
    pub hgvs: bool,
    pub hgvsg: bool,
    /// `--hgvsp_use_prediction`: HGVSp's description after `p.` is parenthesised.
    pub hgvsp_use_prediction: bool,
    /// `--remove_hgvsp_version`: HGVSp's reference name carries no version.
    pub remove_hgvsp_version: bool,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    pub numbers: bool,
    pub domains: bool,
    pub regulatory: bool,
    pub variant_class: bool,
    pub overlaps: bool,
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomade: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub pubmed: bool,
    pub old_maf: bool,
    pub gene_phenotype: bool,
    pub mirna: bool,
    pub check_existing: bool,
    /// `--no_check_alleles`: a known variant co-locates on its span alone.
    pub no_check_alleles: bool,
    /// `--exclude_null_alleles`: a known variant without alleles never co-locates.
    pub exclude_null_alleles: bool,
    /// `--failed <int>`: the highest `failed` value a known variant may carry and
    /// still co-locate; Ensembl VEP's default is 0.
    pub failed: i64,
    /// `--clin_sig_allele`: CLIN_SIG reports the row's allele's own ClinVar
    /// entries (true, Ensembl VEP's default) or every variant's terms (false).
    pub clin_sig_allele: bool,
    pub clinvar_somatic_classification: bool,
    pub var_synonyms: bool,
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
    /// The `--transcript_filter` strings, in command-line order; empty when
    /// none was passed.
    pub transcript_filter: Vec<String>,
    pub allele_number: bool,
    pub show_ref_allele: bool,
    pub minimal: bool,
    pub total_length: bool,

    pub no_stats: bool,
    pub stats_file: Option<String>,
    pub stats_text: bool,

    pub verbose: bool,
    pub quiet: bool,
    /// The `INFO:` lines the resolution owes the run, in Ensembl VEP's order
    /// (`Runner.pm:117-120` before `post_setup_checks`), printed by the runner
    /// unless `--quiet`.
    pub notices: Vec<&'static str>,

    pub plugin: Vec<String>,
    pub custom: Vec<String>,

    pub refseq: bool,
    pub merged: bool,
    pub gencode_basic: bool,
    pub gencode_primary: bool,
    pub flag_gencode_primary: bool,

    pub vcf_info_field: String,
    /// Keep an input VCF's existing annotation beside the new one (`--keep_csq`).
    pub keep_csq: bool,

    pub json_cache: Option<String>,

    /// The `--fields` list: the columns of the tab output and the subfields of
    /// the VCF entry, replacing the whole default list; `None` keeps the defaults.
    pub fields: Option<Vec<String>>,
    pub dont_skip: bool,
    pub allow_non_variant: bool,
}

/// The output format as the writers name it: `OutputFactory.pm:229-231` lower-cases
/// the value and its `%FORMAT_MAP` (`:97-103`) reads `ensembl` as the default `vep`
/// format.
pub fn canonical_output_format(value: &str) -> String {
    match value.to_ascii_lowercase().as_str() {
        "ensembl" => "vep".to_string(),
        lower => lower.to_string(),
    }
}

/// A `@LIST_FLAGS` value split as `Config.pm:757-762` splits it, with Perl's
/// `split`: an inner empty field stays, trailing empty fields are dropped.
fn split_list_flag(value: &str) -> Vec<String> {
    let mut parts: Vec<String> = value.split(',').map(str::to_string).collect();
    while parts.last().is_some_and(String::is_empty) {
        parts.pop();
    }
    parts
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

/// The GENCODE primary set exists for human GRCh38 alone, so `--gencode_primary`
/// and `--flag_gencode_primary` are each warned about when `--assembly` names
/// another assembly (Config.pm:722-729), `hg38` read as the GRCh38 it stands
/// for here. Ensembl VEP reads the assembly before its option sets apply, so
/// `--grch37` alone, which sets the assembly through one, raises no warning; it
/// also prints the text on standard output, where vep-rs routes it through the
/// warning sink.
fn gencode_primary_assembly_warnings(state: &FlagState) -> Vec<String> {
    let other_assembly = state
        .scalar("assembly")
        .map(crate::cache_resolve::canonical_assembly)
        .is_some_and(|assembly| !assembly.eq_ignore_ascii_case("grch38"));
    if !other_assembly {
        return Vec::new();
    }
    ["gencode_primary", "flag_gencode_primary"]
        .into_iter()
        .filter(|flag| state.is_active(flag))
        .map(|flag| {
            format!("--{flag} option is currently only available for human on the GRCh38 assembly")
        })
        .collect()
}

/// `--var_synonyms` switches `check_existing` on when it is off, with the
/// notice Ensembl VEP's runner prints as it does so (`Runner.pm:117-120` at
/// release/116.2, a status message `--quiet` withholds); the notice is returned
/// when the flag was switched on here, `None` otherwise, for [`Config::notices`].
fn var_synonyms_implies_check_existing(state: &mut FlagState) -> Option<&'static str> {
    if !state.is_active("var_synonyms") || state.is_active("check_existing") {
        return None;
    }
    state.set("check_existing", "1");
    Some("INFO: Colocated variants will be provided when using --var_synonyms")
}

impl Config {
    /// Parses a command line (the binary name first) and resolves it.
    pub fn from_argv<I, T>(argv: I) -> Result<Self, ParseFailure>
    where
        I: IntoIterator<Item = T>,
        T: Into<OsString> + Clone,
    {
        Config::from_parsed(parse_command_line(argv)?).map_err(ParseFailure::Refused)
    }

    /// Resolves a parsed command line: the error table's refusals first, then
    /// Ensembl VEP's option sets and option rules on the flag state
    /// (`Config.pm:766-769`), then vep-rs's own rules, then the configuration.
    pub fn from_parsed(parsed: ParsedArgs) -> Result<Self, CliError> {
        let ParsedArgs { args, mut state } = parsed;
        for entry in ERROR_TABLE.iter().filter(|e| e.declared) {
            let Treatment::Refused { component, message } = entry.treatment else {
                continue;
            };
            if state.is_given(entry.flag) {
                return Err(CliError::new(component, message));
            }
        }
        for warning in gencode_primary_assembly_warnings(&state) {
            tracing::warn!("{}", warning);
        }
        if let Some(format) = state.scalar("output_format") {
            let canonical = canonical_output_format(format);
            state.set("output_format", &canonical);
        }
        state.apply_option_sets();
        state
            .check(state.is_active("safe"))
            .map_err(CliError::cli)?;
        for rule in OWN_RULES {
            if (rule.applies)(&state) {
                return Err(CliError::new(rule.component, rule.message));
            }
        }
        let distance =
            parse_distance(&args.distance).map_err(|message| CliError::new("config", message))?;
        // After every refusal: Ensembl VEP's `Runner.pm:117-120` runs once
        // `Config.pm:769` has checked the command line, so a refused one prints
        // its error alone.
        let notices = var_synonyms_implies_check_existing(&mut state)
            .into_iter()
            .collect();
        let mut config = Config::resolve(args, &state, distance);
        config.notices = notices;
        Ok(config)
    }

    /// The configuration from the arguments and the flag state after the option
    /// sets: a value an option set assigns comes from the state, everything else
    /// from the arguments.
    fn resolve(args: Args, state: &FlagState, distance: (u64, u64)) -> Self {
        let on = |flag: &str| state.is_active(flag);
        let output_format = state.scalar("output_format").unwrap_or("vep").to_string();
        // The VCF CSQ layout names SYMBOL, BIOTYPE, EXON and INTRON whatever the
        // flags, so a VCF output fills them by either route; Ensembl VEP's `vcf`
        // option set does so for `--vcf` alone.
        let vcf_output = output_format == "vcf";

        // Getopt::Long accumulates `-o` (`output_file|o=s@`) and the vep script keeps
        // the last (vep:33); `STDOUT` is matched in any case (BaseRunner.pm:213).
        let last = args
            .output_file
            .last()
            .cloned()
            .unwrap_or_else(|| "variant_effect_output.txt".to_string());
        let warning_file = warning_file_for(&args);
        let output_file = if last.eq_ignore_ascii_case("STDOUT") {
            "STDOUT".to_string()
        } else {
            last
        };

        let dir = expand_tilde(&args.dir);
        // A directory given to --cache stands for --dir_cache when that flag is
        // absent (Config.pm:716-718); any other value of --cache only switches it on.
        // The bare flag's stand-in value is never a directory: Perl's bare `cache:s`
        // is the empty string.
        let cache_dir_value = args
            .cache
            .as_deref()
            .filter(|value| *value != CACHE_FLAG_ALONE)
            .map(expand_tilde)
            .filter(|value| std::path::Path::new(value).is_dir());
        let dir_cache = args
            .dir_cache
            .map(|d| expand_tilde(&d))
            .or(cache_dir_value)
            .unwrap_or_else(|| dir.clone());
        let dir_plugins = args
            .dir_plugins
            .map(|d| expand_tilde(&d))
            .unwrap_or_else(|| format!("{dir}/Plugins"));

        let cache_version = args.cache_version.unwrap_or(vep_core::VEP_VERSION);

        let fields = args.fields.as_deref().map(split_list_flag);

        if let Some(notice) = regulatory_unimplemented_notice(args.regulatory, args.everything) {
            tracing::warn!("{}", notice);
        }

        // `shift_hgvs` is a setting only as 0 or 1 (BaseRunner.pm:498), and 0 sets
        // `shift_3prime` to 0 whatever the command line said (Runner.pm:792-795).
        let shift_hgvs = args.shift_hgvs.as_deref() != Some("0");
        let shift_3prime = shift_hgvs && args.shift_3prime.as_deref() == Some("1");

        Config {
            input_file: args.input_file,
            input_data: args.input_data,
            format: state.scalar("format").unwrap_or("guess").to_string(),
            output_file,
            output_format,
            parquet_shape: args.parquet_shape,
            keep_intermediate: args.keep_intermediate,
            parquet_row_group_size: args.parquet_row_group_size,
            compress_output: args.compress_output,
            force_overwrite: args.force_overwrite,
            no_headers: args.no_headers,
            mask_header_cache_path: args.mask_header_cache_path,
            warning_file,
            skipped_variants_file: args.skipped_variants_file,
            cache: on("cache"),
            offline: on("offline"),
            dir,
            dir_cache,
            dir_plugins,
            safe: on("safe"),
            cache_version,
            show_cache_info: args.show_cache_info,
            fasta: args.fasta,
            fasta_dir: args.fasta_dir,
            no_fasta: args.no_fasta,
            synonyms: args.synonyms.as_deref().map(expand_tilde),
            shift_3prime,
            shift_hgvs,
            shift_genomic: args.shift_genomic.as_deref() == Some("1"),
            shift_length: args.shift_length,
            species: args.species,
            assembly: state.scalar("assembly").map(str::to_string),
            buffer_size: args.buffer_size,
            max_sv_size: u64::try_from(args.max_sv_size).ok(),
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
            everything: on("everything"),
            symbol: on("symbol") || vcf_output,
            biotype: on("biotype") || vcf_output,
            canonical: on("canonical"),
            mane: on("mane"),
            mane_select: on("mane_select"),
            tsl: on("tsl"),
            appris: on("appris"),
            ccds: on("ccds"),
            protein: on("protein"),
            protein_version: args.protein_version,
            uniprot: on("uniprot"),
            xref_refseq: args.xref_refseq,
            hgvs: on("hgvs"),
            hgvsg: args.hgvsg,
            hgvsp_use_prediction: args.hgvsp_use_prediction,
            remove_hgvsp_version: args.remove_hgvsp_version,
            sift: state.scalar("sift").map(str::to_string),
            polyphen: state.scalar("polyphen").map(str::to_string),
            numbers: on("numbers") || vcf_output,
            domains: on("domains"),
            regulatory: on("regulatory"),
            variant_class: on("variant_class"),
            overlaps: args.overlaps,
            af: on("af"),
            af_1kg: on("af_1kg"),
            af_gnomade: on("af_gnomade"),
            af_gnomadg: on("af_gnomadg"),
            max_af: on("max_af"),
            pubmed: on("pubmed"),
            old_maf: on("old_maf"),
            gene_phenotype: on("gene_phenotype"),
            mirna: on("mirna"),
            check_existing: on("check_existing"),
            no_check_alleles: args.no_check_alleles,
            exclude_null_alleles: args.exclude_null_alleles,
            failed: args.failed,
            clin_sig_allele: args.clin_sig_allele != 0,
            clinvar_somatic_classification: args.clinvar_somatic_classification,
            var_synonyms: on("var_synonyms"),
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
            pick_order: args.pick_order,
            transcript_filter: args.transcript_filter,
            allele_number: on("allele_number"),
            show_ref_allele: args.show_ref_allele,
            minimal: args.minimal,
            total_length: args.total_length,
            no_stats: on("no_stats"),
            stats_file: args.stats_file,
            stats_text: args.stats_text,
            verbose: args.verbose,
            quiet: args.quiet,
            notices: Vec::new(),
            plugin: args.plugin,
            custom: state.list("custom"),
            refseq: args.refseq,
            merged: args.merged,
            gencode_basic: args.gencode_basic,
            gencode_primary: args.gencode_primary,
            flag_gencode_primary: args.flag_gencode_primary,
            vcf_info_field: args.vcf_info_field,
            keep_csq: args.keep_csq,
            json_cache: args.json_cache.as_deref().map(expand_tilde),
            fields,
            dont_skip: args.dont_skip,
            allow_non_variant: on("allow_non_variant"),
        }
    }
}

/// Auto-detect cap for `--decompression_threads=0`: decompression throughput
/// plateaus between 4 and 10 threads and is flat past that, so the auto default
/// stops at 10 rather than oversubscribing large hosts. User-supplied values are
/// not capped; any positive integer pins the worker count.
const DECOMPRESSION_THREADS_AUTO_CAP: usize = 10;

/// Auto-detect cap for `--fork=0`. One reader thread parses and hands out the
/// batches, so the wall time stops falling past about 32 workers and rises on
/// wider hosts, where the extra workers wait at the batch barrier and contend
/// for the allocator: gnomAD v4.1 chromosome 2 took 1.8 times longer at 192
/// workers than at 32 on a 192-vCPU host. User-supplied values are not capped.
const FORK_AUTO_CAP: usize = 32;

/// Resolve the --fork value. `0` means "auto": use the number of available
/// logical CPUs, capped at [`FORK_AUTO_CAP`]. Any non-zero value is returned
/// as-is so users can still pin the worker count explicitly.
///
/// Public so an embedding consumer of this crate can apply the same
/// auto-detect default without going through the full `Config::from_argv`
/// pipeline.
pub fn resolve_fork(raw: usize) -> usize {
    if raw != 0 {
        return raw;
    }
    std::thread::available_parallelism()
        .map(|n| n.get().min(FORK_AUTO_CAP))
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

/// The `--distance` value as (upstream, downstream) bp: one integer for both, or
/// the pair `up,down`. Ensembl VEP reads the pair with `$d ||= $u` (`BaseRunner.pm`
/// 506), which fires on the false string `0` alone, so the downstream part `0` is
/// the upstream value while `00` or `+0` is a downstream distance of 0. Ensembl VEP
/// numifies a non-numeric part to 0 under an `Argument "x" isn't numeric` warning
/// and drops the up- or downstream rows that part governs (both for a single
/// value); here it is a refusal naming the value, as is a part padded with
/// whitespace, which Perl numifies but `||=` reads as set.
fn parse_distance(s: &str) -> Result<(u64, u64), String> {
    let refused = || format!("--distance expects one or two integers (up[,down]); got {s:?}");
    let integer = |part: &str| part.parse::<u64>().map_err(|_| refused());
    let (up, down) = match s.split_once(',') {
        None => (s, s),
        Some((up, "0")) => (up, up),
        Some((up, down)) => (up, down),
    };
    Ok((integer(up)?, integer(down)?))
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
    use crate::option_rules::OPTION_SETS;
    use clap::CommandFactory;

    fn config(argv: &[&str]) -> Config {
        Config::from_argv(argv).unwrap()
    }

    fn refusal(argv: &[&str]) -> String {
        match Config::from_argv(argv).unwrap_err() {
            ParseFailure::Refused(e) => e.to_string(),
            ParseFailure::Usage(e) => panic!("usage error, not a refusal: {e}"),
        }
    }

    /// A clap usage error's message, for a command line the parser refuses.
    fn usage_error(argv: &[&str]) -> String {
        match Config::from_argv(argv).unwrap_err() {
            ParseFailure::Usage(e) => e.to_string(),
            ParseFailure::Refused(e) => panic!("refused, not a usage error: {e}"),
        }
    }

    /// The auto thread count never exceeds the cap, whatever the host offers; an
    /// explicit value passes through uncapped.
    #[test]
    fn auto_fork_is_capped_and_explicit_fork_is_not() {
        let auto = resolve_fork(0);
        assert!(
            (1..=FORK_AUTO_CAP).contains(&auto),
            "auto fork {auto} outside 1..={FORK_AUTO_CAP}"
        );
        assert_eq!(resolve_fork(1), 1);
        assert_eq!(resolve_fork(64), 64);
    }

    /// The flags `--everything` turns on are Ensembl VEP's set (Config.pm:351-380)
    /// closed under its own option sets (`hgvs` at :382-388, `mane` at :580-586,
    /// `check_existing` at :448-453), and nothing else: `allele_number` and
    /// `hgvsg` stay off.
    #[test]
    fn everything_equals_ensembl_vep_set() {
        let perl_set = [
            "sift",
            "polyphen",
            "ccds",
            "hgvs",
            "symbol",
            "numbers",
            "domains",
            "regulatory",
            "canonical",
            "protein",
            "biotype",
            "af",
            "af_1kg",
            "af_gnomade",
            "af_gnomadg",
            "max_af",
            "pubmed",
            "uniprot",
            "mane",
            "tsl",
            "appris",
            "variant_class",
            "gene_phenotype",
            "mirna",
        ];
        let everything = OPTION_SETS
            .iter()
            .find(|s| s.flags == ["everything"])
            .unwrap();
        let keys: Vec<&str> = everything.set.iter().map(|(k, _)| *k).collect();
        assert_eq!(keys, perl_set);
        let mut implied = FlagState::implied_by("everything");
        implied.sort_unstable();
        let mut expected: Vec<&str> = perl_set.to_vec();
        expected.extend([
            "hgvsc",
            "hgvsp",
            "mane_select",
            "mane_plus_clinical",
            "check_existing",
        ]);
        expected.sort_unstable();
        assert_eq!(implied, expected);

        let c = config(&["vep", "--everything"]);
        let flags: [(&str, bool); 29] = [
            ("ccds", c.ccds),
            ("hgvs", c.hgvs),
            ("symbol", c.symbol),
            ("numbers", c.numbers),
            ("domains", c.domains),
            ("regulatory", c.regulatory),
            ("canonical", c.canonical),
            ("protein", c.protein),
            ("biotype", c.biotype),
            ("af", c.af),
            ("af_1kg", c.af_1kg),
            ("af_gnomade", c.af_gnomade),
            ("af_gnomadg", c.af_gnomadg),
            ("max_af", c.max_af),
            ("pubmed", c.pubmed),
            ("uniprot", c.uniprot),
            ("mane", c.mane),
            ("mane_select", c.mane_select),
            ("tsl", c.tsl),
            ("appris", c.appris),
            ("variant_class", c.variant_class),
            ("gene_phenotype", c.gene_phenotype),
            ("mirna", c.mirna),
            ("check_existing", c.check_existing),
            ("hgvsg", c.hgvsg),
            ("allele_number", c.allele_number),
            ("xref_refseq", c.xref_refseq),
            ("protein_version", c.protein_version),
            ("minimal", c.minimal),
        ];
        for (flag, on) in flags {
            assert_eq!(on, implied.contains(&flag), "{flag} under --everything");
        }
        assert_eq!(c.sift.as_deref(), Some("b"));
        assert_eq!(c.polyphen.as_deref(), Some("b"));
        let plain = config(&["vep"]);
        assert!(!plain.symbol && !plain.hgvs && !plain.mane_select && plain.sift.is_none());
    }

    /// The `--everything` help, which the CLI reference prints, states the one
    /// exclusion: the regulatory columns are not produced.
    #[test]
    fn everything_documents_the_regulatory_exclusion() {
        let mut command = Args::command_with_error_table();
        command.build();
        let help = command
            .get_arguments()
            .find(|a| a.get_long() == Some("everything"))
            .and_then(|a| a.get_long_help().or(a.get_help()))
            .map(|h| h.to_string())
            .unwrap();
        assert!(
            help.contains("except the regulatory columns")
                && help.contains("no RegulatoryFeature or MotifFeature row is written"),
            "{help}"
        );
    }

    /// An option set overwrites an explicit value, as `apply_option_sets` does
    /// (Config.pm:875): `--everything --sift p` runs SIFT in `b` mode.
    #[test]
    fn option_sets_override_explicit_values_as_perl_does() {
        let c = config(&["vep", "--everything", "--sift", "p"]);
        assert_eq!(c.sift.as_deref(), Some("b"));
        let c = config(&["vep", "--sift", "p"]);
        assert_eq!(c.sift.as_deref(), Some("p"));
    }

    #[test]
    fn test_offline_implies_cache() {
        let c = config(&["vep", "--offline"]);
        assert!(c.cache);
        assert!(c.offline);
    }

    /// `--vcf` switches on symbol, biotype and numbers (the option set at
    /// Config.pm:500-508); `--output_format vcf` does the same here, because the
    /// CSQ layout names SYMBOL, BIOTYPE, EXON and INTRON whatever the flags.
    #[test]
    fn test_vcf_shorthand() {
        for argv in [
            &["vep", "--vcf"][..],
            &["vep", "--output_format", "vcf"][..],
        ] {
            let c = config(argv);
            assert_eq!(c.output_format, "vcf", "{argv:?}");
            assert!(c.symbol && c.biotype && c.numbers, "{argv:?}");
        }
        let plain = config(&["vep"]);
        assert!(!plain.symbol && !plain.biotype && !plain.numbers);
    }

    #[test]
    fn test_json_shorthand() {
        assert_eq!(config(&["vep", "--json"]).output_format, "json");
    }

    #[test]
    fn test_tab_shorthand() {
        assert_eq!(config(&["vep", "--tab"]).output_format, "tab");
    }

    /// `OutputFactory.pm:229-231` lower-cases `--output_format` and its `%FORMAT_MAP`
    /// (`:97-103`) names `ensembl` and `vep` as one format; a value outside the map is
    /// refused (`Unknown or unsupported output format`), which here is the value
    /// enumeration's syntax error listing the six values, `parquet` being vep-rs's
    /// own. The canonical name is what the option rules read, so `PARQUET` meets the
    /// Parquet rules.
    #[test]
    fn output_format_is_lower_cased_validated_and_ensembl_is_vep() {
        assert_eq!(
            config(&["vep", "--output_format", "VCF"]).output_format,
            "vcf"
        );
        assert_eq!(
            config(&["vep", "--output_format", "Json"]).output_format,
            "json"
        );
        for spelling in ["ensembl", "Ensembl", "ENSEMBL", "vep", "VEP"] {
            assert_eq!(
                config(&["vep", "--output_format", spelling]).output_format,
                "vep",
                "{spelling}"
            );
        }
        let message = usage_error(&["vep", "--output_format", "xml"]);
        assert!(
            message.contains("invalid value 'xml' for '--output_format")
                && message.contains("[possible values: vep, ensembl, vcf, json, tab, parquet]"),
            "{message}"
        );
        assert_eq!(
            refusal(&["vep", "--output_format", "PARQUET", "-o", "STDOUT"]),
            "ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT"
        );
        // `--output_format VCF` fills the CSQ columns as `--vcf` does.
        assert!(config(&["vep", "--output_format", "VCF"]).symbol);
    }

    /// `--json` refuses `--vcf` and `--tab` (Config.pm:623); the pair is named in the
    /// incompatibility table's order, `json` before `vcf` and `tab`, where Ensembl
    /// VEP iterates a hash (`:962`) and names either flag first from run to run.
    #[test]
    fn json_with_vcf_or_tab_is_perls_error() {
        assert_eq!(
            refusal(&["vep", "--json", "--vcf"]),
            "ERROR: [cli] Can't use --json and --vcf together"
        );
        assert_eq!(
            refusal(&["vep", "--json", "--tab"]),
            "ERROR: [cli] Can't use --json and --tab together"
        );
        assert_eq!(config(&["vep", "--json"]).output_format, "json");
    }

    /// `--vcf` refuses `--json` and `--tab` (Config.pm:624). With `--json` the `json`
    /// row of the table fires first, so the line reads `--json and --vcf`.
    #[test]
    fn vcf_with_json_or_tab_is_perls_error() {
        assert_eq!(
            refusal(&["vep", "--vcf", "--json"]),
            "ERROR: [cli] Can't use --json and --vcf together"
        );
        assert_eq!(
            refusal(&["vep", "--vcf", "--tab"]),
            "ERROR: [cli] Can't use --vcf and --tab together"
        );
    }

    /// `--most_severe` lists `vcf` among its incompatibilities (Config.pm:616), and
    /// the `vcf` option set has switched on `biotype`, the row's first entry, by the
    /// time the checks run (`:766-769`), so the line names that pair; `--most_severe`
    /// alone with `--vcf` cannot be reached. The rule keys on the flag: Ensembl VEP
    /// accepts `--output_format vcf --most_severe`, and so does vep-rs.
    #[test]
    fn vcf_with_most_severe_is_perls_error() {
        assert_eq!(
            refusal(&["vep", "--vcf", "--most_severe"]),
            "ERROR: [cli] Can't use --most_severe and --biotype together"
        );
        let mut state = FlagState::default();
        state.set("most_severe", "1");
        state.set("vcf", "1");
        assert_eq!(
            state.check(false).unwrap_err(),
            "Can't use --most_severe and --vcf together"
        );
        assert!(Config::from_argv(["vep", "--output_format", "vcf", "--most_severe"]).is_ok());
    }

    /// `--tab` refuses `--vcf` and `--json` (Config.pm:625); both other rows precede
    /// the `tab` row, so each line names the other flag first.
    #[test]
    fn tab_with_vcf_or_json_is_perls_error() {
        assert_eq!(
            refusal(&["vep", "--tab", "--vcf"]),
            "ERROR: [cli] Can't use --vcf and --tab together"
        );
        assert_eq!(
            refusal(&["vep", "--tab", "--json"]),
            "ERROR: [cli] Can't use --json and --tab together"
        );
    }

    /// `--fields` has no Parquet column set to replace: the pair is refused under the
    /// output component, while `--tab` and `--vcf` take the list.
    #[test]
    fn parquet_with_fields_is_refused() {
        assert_eq!(
            refusal(&[
                "vep",
                "--output_format",
                "parquet",
                "--fields",
                "Allele,Gene",
                "-o",
                "x.parquet"
            ]),
            "ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed"
        );
        assert_eq!(
            config(&["vep", "--tab", "--fields", "Allele,Gene"]).fields,
            Some(vec!["Allele".to_string(), "Gene".to_string()])
        );
        assert!(Config::from_argv(["vep", "--vcf", "--fields", "Allele"]).is_ok());
    }

    /// Parquet output is compressed by its own codec, so `--compress_output` is
    /// refused with it and accepted with every text format.
    #[test]
    fn parquet_with_compress_output_is_refused() {
        assert_eq!(
            refusal(&[
                "vep",
                "--output_format",
                "parquet",
                "--compress_output",
                "gzip",
                "-o",
                "x.parquet"
            ]),
            "ERROR: [output] --compress_output applies to the default, tab, VCF and JSON outputs; Parquet output is compressed by its own codec"
        );
        for format in ["vep", "tab", "vcf", "json"] {
            let c = config(&[
                "vep",
                "--output_format",
                format,
                "--compress_output",
                "bgzip",
            ]);
            assert_eq!(c.compress_output.as_deref(), Some("bgzip"), "{format}");
        }
    }

    /// `%VALID` holds `gzip` and `bgzip` for `--compress_output` (Config.pm:603) and
    /// `check_config` (`:945-957`) refuses another value with Perl's two-line
    /// message; the declared flag carries the same list as a value enumeration, so
    /// on the command line the refusal is the syntax error listing both values.
    #[test]
    fn compress_output_rejects_unknown_value_with_perls_message() {
        let mut state = FlagState::default();
        state.set("compress_output", "zip");
        assert_eq!(
            state.check(false).unwrap_err(),
            "\"zip\" is not a valid value for --compress_output\nValid values: gzip, bgzip"
        );
        let message = usage_error(&["vep", "--compress_output", "zip"]);
        assert!(
            message.contains("invalid value 'zip' for '--compress_output")
                && message.contains("[possible values: gzip, bgzip]"),
            "{message}"
        );
        assert_eq!(
            config(&["vep", "--compress_output", "gzip"])
                .compress_output
                .as_deref(),
            Some("gzip")
        );
        assert_eq!(config(&["vep"]).compress_output, None);
    }

    /// `--fields` is a `@LIST_FLAGS` value (Config.pm:327-336), split on commas as
    /// Perl's `split` does: nothing is trimmed, an inner empty field stays and a
    /// trailing comma opens no field.
    #[test]
    fn fields_split_follows_perl() {
        let fields = |value: &str| config(&["vep", "--tab", "--fields", value]).fields.unwrap();
        assert_eq!(fields("Allele,Gene,"), ["Allele", "Gene"]);
        assert_eq!(fields("Allele, Gene"), ["Allele", " Gene"]);
        assert_eq!(fields("Allele,,Gene"), ["Allele", "", "Gene"]);
        assert_eq!(config(&["vep", "--tab"]).fields, None);
    }

    /// `--keep_csq` is read; `--keep_ann`, its synonym Ensembl VEP never reads, is not.
    #[test]
    fn keep_csq_is_read_and_keep_ann_is_not() {
        assert!(config(&["vep", "--vcf", "--keep_csq"]).keep_csq);
        assert!(!config(&["vep", "--vcf", "--keep_ann"]).keep_csq);
        assert!(!config(&["vep", "--vcf"]).keep_csq);
    }

    /// `--minimal` switches on `--allele_number` (the option set at Config.pm:524-529).
    #[test]
    fn minimal_implies_allele_number() {
        assert!(config(&["vep", "--minimal"]).allele_number);
        assert!(!config(&["vep"]).allele_number);
    }

    /// `--var_synonyms` without `--check_existing` switches the latter on and owes
    /// the notice `Runner.pm:117-120` (release/116.2) prints as it does so; given
    /// beside `--check_existing`, or with the flag off, nothing is switched and
    /// no notice is due. The notice is a status message, so `--quiet` withholds
    /// it (BaseVEP.pm:782) while the implication stands. A refused command line
    /// is refused as it is without the flag: `Config.pm:769` checks it before
    /// the runner's implication runs.
    #[test]
    fn var_synonyms_implies_check_existing_with_info_line() {
        const NOTICE: &str = "INFO: Colocated variants will be provided when using --var_synonyms";
        let implied = config(&["vep", "--var_synonyms"]);
        assert!(implied.var_synonyms && implied.check_existing);
        assert_eq!(implied.notices, [NOTICE]);
        let quiet = config(&["vep", "--var_synonyms", "--quiet"]);
        assert!(quiet.var_synonyms && quiet.check_existing);
        assert_eq!(quiet.notices, [NOTICE], "the runner withholds it");
        let off = config(&["vep", "--check_existing"]);
        assert!(!off.var_synonyms && off.check_existing);
        assert!(off.notices.is_empty());
        assert!(!config(&["vep"]).check_existing);
        assert!(config(&["vep", "--var_synonyms", "--check_existing"])
            .notices
            .is_empty());
        assert_eq!(
            refusal(&["vep", "--var_synonyms", "--json", "--vcf"]),
            "ERROR: [cli] Can't use --json and --vcf together"
        );

        let state_after_option_sets = |argv: &[&str]| {
            let mut state = parse_command_line(argv).unwrap().state;
            state.apply_option_sets();
            state
        };
        let mut state = state_after_option_sets(&["vep", "--var_synonyms"]);
        assert_eq!(
            var_synonyms_implies_check_existing(&mut state),
            Some(NOTICE)
        );
        assert!(state.is_active("check_existing"));
        let mut state = state_after_option_sets(&["vep", "--var_synonyms", "--check_existing"]);
        assert_eq!(var_synonyms_implies_check_existing(&mut state), None);
        let mut state = state_after_option_sets(&["vep", "--pubmed"]);
        assert_eq!(var_synonyms_implies_check_existing(&mut state), None);
    }

    /// `--clin_sig_allele` is Ensembl VEP's `clin_sig_allele=i` (Config.pm:151),
    /// on by default (:316): any integer parses and `0` alone switches the
    /// allele-specific significance off; `--clinvar_somatic_classification` and
    /// `--var_synonyms` are bare flags (:152, :230).
    #[test]
    fn clin_sig_allele_takes_an_integer_and_defaults_on() {
        assert!(config(&["vep"]).clin_sig_allele);
        assert!(!config(&["vep", "--clin_sig_allele", "0"]).clin_sig_allele);
        assert!(config(&["vep", "--clin_sig_allele", "1"]).clin_sig_allele);
        assert!(config(&["vep", "--clin_sig_allele=2"]).clin_sig_allele);
        assert!(usage_error(&["vep", "--clin_sig_allele", "yes"]).contains("invalid value"));
        let flags = config(&["vep", "--clinvar_somatic_classification", "--var_synonyms"]);
        assert!(flags.clinvar_somatic_classification && flags.var_synonyms);
        assert!(!config(&["vep"]).clinvar_somatic_classification);
    }

    /// Two output-format shorthands together are Ensembl VEP's incompatibility
    /// (Config.pm:623-625), as are `--quiet --verbose` (:620) and `--refseq --merged`
    /// (:621).
    #[test]
    fn incompatible_pairs_are_refused() {
        assert_eq!(
            refusal(&["vep", "--vcf", "--json"]),
            "ERROR: [cli] Can't use --json and --vcf together"
        );
        assert_eq!(
            refusal(&["vep", "--quiet", "--verbose"]),
            "ERROR: [cli] Can't use --quiet and --verbose together"
        );
        assert_eq!(
            refusal(&["vep", "--refseq", "--merged"]),
            "ERROR: [cli] Can't use --refseq and --merged together"
        );
    }

    /// `--most_severe` refuses every flag of its `%INCOMPATIBLE` list
    /// (Config.pm:616) with Ensembl VEP's message; `--per_gene`, which the list
    /// omits, is accepted beside it.
    #[test]
    fn most_severe_refuses_its_incompatible_flags() {
        let incompatible = [
            "biotype",
            "no_intergenic",
            "protein",
            "symbol",
            "coding_only",
            "ccds",
            "mane",
            "mane_select",
            "canonical",
            "xref_refseq",
            "numbers",
            "domains",
            "tsl",
            "appris",
            "uniprot",
            "pick",
            "flag_pick",
            "pick_allele",
            "flag_pick_allele",
            "summary",
        ];
        for other in incompatible {
            let flag = format!("--{other}");
            assert_eq!(
                refusal(&["vep", "--most_severe", &flag]),
                format!("ERROR: [cli] Can't use --most_severe and --{other} together"),
                "{flag}"
            );
        }
        // The `vcf` option set switches on `biotype` (Config.pm:500-508), the
        // first flag of the list, so that pair is named.
        assert_eq!(
            refusal(&["vep", "--most_severe", "--vcf"]),
            "ERROR: [cli] Can't use --most_severe and --biotype together"
        );
        for (tool, value) in [("sift", "b"), ("polyphen", "p")] {
            let flag = format!("--{tool}");
            assert_eq!(
                refusal(&["vep", "--most_severe", &flag, value]),
                format!("ERROR: [cli] Can't use --most_severe and --{tool} together")
            );
        }
        let accepted = config(&["vep", "--most_severe", "--per_gene"]);
        assert!(accepted.most_severe && accepted.per_gene);
    }

    /// `--summary` refuses the same list with `--most_severe` in place of
    /// `--summary` (Config.pm:617); the pair of the two is named by the
    /// `most_severe` row, which the table lists first.
    #[test]
    fn summary_refuses_its_incompatible_flags() {
        for other in ["biotype", "coding_only", "pick", "flag_pick_allele"] {
            let flag = format!("--{other}");
            assert_eq!(
                refusal(&["vep", "--summary", &flag]),
                format!("ERROR: [cli] Can't use --summary and --{other} together"),
                "{flag}"
            );
        }
        assert_eq!(
            refusal(&["vep", "--summary", "--vcf"]),
            "ERROR: [cli] Can't use --summary and --biotype together"
        );
        assert_eq!(
            refusal(&["vep", "--summary", "--most_severe"]),
            "ERROR: [cli] Can't use --most_severe and --summary together"
        );
        let accepted = config(&["vep", "--summary", "--pick_allele_gene"]);
        assert!(accepted.summary && accepted.pick_allele_gene);
    }

    /// `--pick_order` takes Ensembl VEP's `%VALID` tokens (Config.pm:601) in any
    /// order and count; a token outside the list is refused when the line is
    /// parsed, with the valid values listed, as every `%VALID` flag is.
    #[test]
    fn pick_order_rejects_a_bad_token_listing_perls_valid_values() {
        assert_eq!(
            config(&["vep", "--pick_order", "rank,canonical,length"])
                .pick_order
                .as_deref(),
            Some(
                &[
                    "rank".to_string(),
                    "canonical".to_string(),
                    "length".to_string()
                ][..]
            )
        );
        assert_eq!(config(&["vep"]).pick_order, None);
        let ParseFailure::Usage(error) =
            Config::from_argv(["vep", "--pick_order", "canonical,x"]).unwrap_err()
        else {
            panic!("a bad token is refused by the parser");
        };
        assert_eq!(error.kind(), clap::error::ErrorKind::InvalidValue);
        let text = error.to_string();
        assert!(text.contains("'x'"), "{text}");
        assert!(
            text.contains(
                "mane_select, mane_plus_clinical, canonical, appris, tsl, biotype, ccds, rank, length, ensembl, refseq"
            ),
            "{text}"
        );
    }

    /// `%INCOMPATIBLE` pairs `gencode_primary` with `gencode_basic` and `refseq`
    /// (Config.pm:621-622), so the two GENCODE sets cannot be combined; either
    /// alone, and `--flag_gencode_primary` beside `--gencode_basic`, resolve.
    #[test]
    fn gencode_primary_refuses_gencode_basic() {
        assert_eq!(
            refusal(&["vep", "--gencode_primary", "--gencode_basic"]),
            "ERROR: [cli] Can't use --gencode_primary and --gencode_basic together"
        );
        assert_eq!(
            refusal(&["vep", "--gencode_primary", "--refseq"]),
            "ERROR: [cli] Can't use --gencode_primary and --refseq together"
        );
        assert!(config(&["vep", "--gencode_basic"]).gencode_basic);
        assert!(config(&["vep", "--gencode_primary"]).gencode_primary);
        let flagged = config(&["vep", "--flag_gencode_primary", "--gencode_basic"]);
        assert!(flagged.flag_gencode_primary && flagged.gencode_basic);
    }

    /// `--flag_gencode_primary` switches the `GENCODE_PRIMARY` column on without
    /// restricting the transcript set, and `--gencode_primary` restricts the set
    /// without adding the column (Config.pm:100, 116: two flags, no option set
    /// joining them).
    #[test]
    fn flag_gencode_primary_sets_column() {
        let column = config(&["vep", "--flag_gencode_primary"]);
        assert!(column.flag_gencode_primary);
        assert!(!column.gencode_primary);
        let set = config(&["vep", "--gencode_primary"]);
        assert!(set.gencode_primary);
        assert!(!set.flag_gencode_primary);
        assert!(FlagState::implied_by("flag_gencode_primary").is_empty());
        assert!(FlagState::implied_by("gencode_primary").is_empty());
    }

    /// Each GENCODE primary flag is warned about when `--assembly` names an
    /// assembly other than GRCh38, in Ensembl VEP's words (Config.pm:722-729),
    /// read before the option sets apply (`:766`), so `--grch37` alone raises
    /// none; the `hg38` and `hg19` spellings count as the assemblies they name.
    #[test]
    fn gencode_primary_warns_off_grch38() {
        let warnings = |argv: &[&str]| {
            gencode_primary_assembly_warnings(&parse_command_line(argv).unwrap().state)
        };
        assert_eq!(
            warnings(&["vep", "--gencode_primary", "--assembly", "GRCh37"]),
            ["--gencode_primary option is currently only available for human on the GRCh38 assembly"]
        );
        assert_eq!(
            warnings(&["vep", "--flag_gencode_primary", "--assembly", "grch37"]),
            ["--flag_gencode_primary option is currently only available for human on the GRCh38 assembly"]
        );
        assert_eq!(
            warnings(&[
                "vep",
                "--gencode_primary",
                "--flag_gencode_primary",
                "--assembly",
                "GRCh37"
            ])
            .len(),
            2
        );
        assert!(warnings(&["vep", "--gencode_primary", "--assembly", "grch38"]).is_empty());
        assert!(warnings(&["vep", "--gencode_primary", "--assembly", "hg38"]).is_empty());
        assert_eq!(
            warnings(&["vep", "--gencode_primary", "--assembly", "hg19"]).len(),
            1
        );
        assert!(warnings(&["vep", "--gencode_primary"]).is_empty());
        assert!(warnings(&["vep", "--gencode_primary", "--grch37"]).is_empty());
        assert!(warnings(&["vep", "--assembly", "GRCh37"]).is_empty());
    }

    /// Every refused entry of the error table is refused with exactly its line,
    /// through the parser for a declared flag and the pre-parse scan otherwise;
    /// every accepted entry parses.
    #[test]
    fn error_table_entries_behave_as_written() {
        let mut command = Args::command();
        command.build();
        let declared: Vec<(String, bool)> = command
            .get_arguments()
            .filter_map(|a| {
                a.get_long().map(|l| {
                    (
                        l.to_string(),
                        matches!(a.get_action(), clap::ArgAction::SetTrue),
                    )
                })
            })
            .collect();
        for entry in ERROR_TABLE {
            let takes_value = declared
                .iter()
                .find(|(l, _)| l == entry.flag)
                .map(|(_, is_flag)| !is_flag);
            assert_eq!(
                takes_value.is_some(),
                entry.declared,
                "--{} declared as the table says",
                entry.flag
            );
            let flag = format!("--{}", entry.flag);
            let mut argv = vec!["vep", flag.as_str()];
            if takes_value == Some(true) {
                argv.push("x");
            }
            match entry.treatment {
                Treatment::Refused { component, message } => {
                    assert_eq!(
                        refusal(&argv),
                        format!("ERROR: [{component}] {message}"),
                        "{flag}"
                    );
                }
                Treatment::AcceptedWithoutEffect { .. } => {
                    assert!(Config::from_argv(&argv).is_ok(), "{flag}");
                }
            }
        }
        assert_eq!(
            refusal(&["vep", "--af_exac=1"]),
            "ERROR: [cli] --af_exac was removed in Ensembl VEP release 115: use --af_gnomade"
        );
    }

    /// `--convert`'s refusal is Ensembl VEP's deprecation text, the same line
    /// whether the table or the option rules produce it.
    #[test]
    fn convert_refusal_is_the_deprecation_text() {
        let entry = error_table_entry("convert").unwrap();
        let Treatment::Refused { message, .. } = entry.treatment else {
            panic!("convert is refused");
        };
        let mut state = FlagState::default();
        state.set("convert", "vcf");
        assert_eq!(state.check(false).unwrap_err(), message);
    }

    /// `custom` is read from the state after the option sets, which the `gff`,
    /// `gtf`, `bigwig`, `phyloP` and `phastCons` sets append to (Config.pm:531-563).
    #[test]
    fn custom_follows_the_state_the_option_sets_append_to() {
        let mut parsed = parse_command_line(["vep", "--custom", "x.bed,X,bed"]).unwrap();
        parsed.state.set("gff", "a.gff");
        let c = Config::from_parsed(parsed).unwrap();
        assert_eq!(c.custom, ["x.bed,X,bed", "file=a.gff,format=gff"]);
        assert_eq!(
            config(&["vep", "--custom", "y.bed,Y,bed"]).custom,
            ["y.bed,Y,bed"]
        );
        assert!(config(&["vep"]).custom.is_empty());
    }

    /// Parquet cannot stream to STDOUT: the output rule refuses the pair.
    #[test]
    fn parquet_to_stdout_is_refused() {
        assert_eq!(
            refusal(&["vep", "--output_format", "parquet", "-o", "stdout"]),
            "ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT"
        );
        assert!(
            Config::from_argv(["vep", "--output_format", "parquet", "-o", "x.parquet"]).is_ok()
        );
    }

    /// `-o` keeps its last value and `STDOUT` is recognised in any case
    /// (BaseRunner.pm:213 `uc eq 'STDOUT'`).
    #[test]
    fn output_file_keeps_the_last_value_and_stdout_in_any_case() {
        assert_eq!(config(&["vep", "-o", "a", "-o", "b"]).output_file, "b");
        for spelling in ["stdout", "StdOut", "STDOUT"] {
            assert_eq!(config(&["vep", "-o", spelling]).output_file, "STDOUT");
        }
        assert_eq!(config(&["vep"]).output_file, "variant_effect_output.txt");
    }

    /// `--grch37` sets the assembly, over an explicit `--assembly`.
    #[test]
    fn grch37_sets_the_assembly() {
        assert_eq!(
            config(&["vep", "--grch37"]).assembly.as_deref(),
            Some("GRCh37")
        );
        assert_eq!(
            config(&["vep", "--assembly", "GRCh38", "--grch37"])
                .assembly
                .as_deref(),
            Some("GRCh37")
        );
        assert_eq!(config(&["vep"]).assembly, None);
    }

    /// A syntax error is clap's, not a refusal.
    #[test]
    fn unknown_flag_is_a_usage_error() {
        match Config::from_argv(["vep", "--no_such_flag"]).unwrap_err() {
            ParseFailure::Usage(e) => {
                assert_eq!(e.kind(), clap::error::ErrorKind::UnknownArgument)
            }
            ParseFailure::Refused(e) => panic!("refused instead: {e}"),
        }
    }

    /// Ensembl VEP splits `--distance` on `,` (`Config.pm` 334) and reads the pair with
    /// `($u, $d) = @$distances; $d ||= $u` (`BaseRunner.pm` 505-506): one value serves
    /// both directions, `5000,2000` is asymmetric, the downstream part `0` is the
    /// upstream value (`2000,0` is `2000,2000`), an upstream 0 stands (`0,2000`), and
    /// `00` or `+0` is a true string in Perl, kept and numified to a downstream 0.
    #[test]
    fn distance_pair_with_zero_downstream_uses_upstream() {
        assert_eq!(parse_distance("5000"), Ok((5000, 5000)));
        assert_eq!(parse_distance("5000,2000"), Ok((5000, 2000)));
        assert_eq!(parse_distance("2000,0"), Ok((2000, 2000)));
        assert_eq!(parse_distance("2000,00"), Ok((2000, 0)));
        assert_eq!(parse_distance("2000,+0"), Ok((2000, 0)));
        assert_eq!(parse_distance("0,2000"), Ok((0, 2000)));
        assert_eq!(parse_distance("0"), Ok((0, 0)));
        assert_eq!(parse_distance("00"), Ok((0, 0)));
        assert_eq!(
            config(&["vep", "--distance", "2000,0"]).distance,
            (2000, 2000)
        );
    }

    /// A part that is not an integer, an empty part, or a third part is refused
    /// with the contract line `ERROR: [config] <message>` naming the value as
    /// written. Ensembl VEP numifies a non-numeric part to 0 under an `Argument "x"
    /// isn't numeric` warning, dropping the up- or downstream rows that part governs,
    /// reads `5000,` as `5000` (Perl's `split` drops a trailing empty field), ignores
    /// a third part, and reads ` 0 ` as a set downstream value (numified to 0) where a
    /// trimmed `0` would read as the upstream value, so a padded part is refused as
    /// well. A negative value in its own token (`--distance -1`) never reaches here:
    /// clap reads `-1` as an unknown flag and exits 2; `--distance=-1` is refused here.
    #[test]
    fn distance_non_numeric_is_an_error() {
        for value in [
            "x", "5000,x", "x,5000", "", "5000,", ",5000", "1.5", "1,2,3", " 2000", "2000, 0",
        ] {
            assert_eq!(
                refusal(&["vep", "--distance", value]),
                format!(
                    "ERROR: [config] --distance expects one or two integers (up[,down]); got {value:?}"
                )
            );
        }
        assert!(matches!(
            Config::from_argv(["vep", "--distance", "-1"]),
            Err(ParseFailure::Usage(_))
        ));
        assert_eq!(
            refusal(&["vep", "--distance=-1"]),
            "ERROR: [config] --distance expects one or two integers (up[,down]); got \"-1\""
        );
    }

    #[test]
    fn test_dir_cache_defaults_to_dir() {
        assert_eq!(
            config(&["vep", "--dir", "/data/vep"]).dir_cache,
            "/data/vep"
        );
    }

    #[test]
    fn test_dir_cache_explicit() {
        assert_eq!(
            config(&["vep", "--dir", "/data/vep", "--dir_cache", "/cache/vep"]).dir_cache,
            "/cache/vep"
        );
    }

    /// A directory given to `--cache` is the cache root when `--dir_cache` is
    /// absent (Config.pm:717); `--dir_cache` wins over it, and a value that is not
    /// a directory leaves the root at `--dir`.
    #[test]
    fn cache_directory_value_stands_for_dir_cache() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().to_str().unwrap();
        assert_eq!(config(&["vep", "--cache", path]).dir_cache, path);
        assert_eq!(
            config(&["vep", "--cache", path, "--dir_cache", "/elsewhere"]).dir_cache,
            "/elsewhere"
        );
        assert_eq!(
            config(&["vep", "--cache", "no-such-directory", "--dir", "/data/vep"]).dir_cache,
            "/data/vep"
        );
        let plain = config(&["vep", "--cache", "--dir", "/data/vep"]);
        assert!(plain.cache);
        assert_eq!(plain.dir_cache, "/data/vep");
        // The bare flag's stand-in value names no directory even when one of that
        // name exists in the working directory.
        struct RemoveOnDrop(Option<std::path::PathBuf>);
        impl Drop for RemoveOnDrop {
            fn drop(&mut self) {
                if let Some(path) = &self.0 {
                    let _ = std::fs::remove_dir(path);
                }
            }
        }
        let marker = std::env::current_dir().unwrap().join(CACHE_FLAG_ALONE);
        // A directory an interrupted test run left behind is reused.
        let created_here = !marker.exists();
        std::fs::create_dir_all(&marker).expect("a directory named for the stand-in value");
        let _marker = RemoveOnDrop(created_here.then_some(marker));
        assert_eq!(
            config(&["vep", "--cache", "--dir", "/data/vep"]).dir_cache,
            "/data/vep"
        );
    }

    /// `--json_cache` is tilde-expanded like the other cache paths.
    #[test]
    fn json_cache_is_tilde_expanded() {
        if let Ok(home) = std::env::var("HOME") {
            assert_eq!(
                config(&["vep", "--json_cache", "~/caches/116"]).json_cache,
                Some(format!("{home}/caches/116"))
            );
        }
        assert_eq!(
            config(&["vep", "--json_cache", "/caches/116"]).json_cache,
            Some("/caches/116".to_string())
        );
    }

    /// Every component an `ERROR: [<component>]` line can name has one row, so
    /// the rendered reference carries one meaning per component.
    #[test]
    fn components_are_unique() {
        let names: std::collections::BTreeSet<&str> =
            COMPONENTS.iter().map(|(name, _)| *name).collect();
        assert_eq!(names.len(), COMPONENTS.len(), "{COMPONENTS:?}");
        for component in ["cli", "output", "cache", "fasta", "filter", "config", "run"] {
            assert!(names.contains(component), "{component}");
        }
    }

    /// `--shift_hgvs 0` sets `shift_3prime` to 0 whatever the command line said
    /// (Runner.pm:792-795); `--shift_hgvs 1` leaves it alone, and a bare run
    /// shifts the HGVS notations and nothing else (BaseRunner.pm:498 reads no
    /// setting from an absent flag).
    #[test]
    fn shift_hgvs_0_forces_shift_3prime_0() {
        let c = config(&["vep", "--shift_hgvs", "0", "--shift_3prime", "0"]);
        assert!(!c.shift_hgvs && !c.shift_3prime);
        let c = config(&["vep", "--shift_hgvs", "0", "--shift_3prime", "1"]);
        assert!(!c.shift_hgvs && !c.shift_3prime);
        let c = config(&["vep", "--shift_hgvs", "0"]);
        assert!(!c.shift_hgvs && !c.shift_3prime);
        let c = config(&["vep", "--shift_hgvs", "1"]);
        assert!(c.shift_hgvs && !c.shift_3prime);
        let c = config(&[
            "vep",
            "--shift_3prime",
            "1",
            "--shift_genomic",
            "1",
            "--shift_length",
        ]);
        assert!(c.shift_hgvs && c.shift_3prime && c.shift_genomic && c.shift_length);
        let plain = config(&["vep"]);
        assert!(plain.shift_hgvs && !plain.shift_3prime && !plain.shift_genomic);
        assert!(!plain.shift_length);
        // The pair Ensembl VEP refuses (Config.pm:630): an active --shift_hgvs
        // beside an active --shift_3prime or --shift_genomic; a 0 is not active.
        assert_eq!(
            refusal(&["vep", "--shift_hgvs", "1", "--shift_3prime", "1"]),
            "ERROR: [cli] Can't use --shift_hgvs and --shift_3prime together"
        );
        assert_eq!(
            refusal(&["vep", "--shift_hgvs", "1", "--shift_genomic", "1"]),
            "ERROR: [cli] Can't use --shift_hgvs and --shift_genomic together"
        );
        assert!(Config::from_argv(["vep", "--shift_hgvs", "0", "--shift_genomic", "0"]).is_ok());
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
    fn regulatory_notice_fires_only_when_requested() {
        assert!(regulatory_unimplemented_notice(false, false).is_none());
        assert!(regulatory_unimplemented_notice(true, false)
            .unwrap()
            .starts_with("--regulatory has no effect"));
        assert!(regulatory_unimplemented_notice(false, true)
            .unwrap()
            .starts_with("--everything implies --regulatory"));
    }

    /// A run error that is not a refusal prints under the `run` component with
    /// its cause chain; a refusal prints as itself.
    #[test]
    fn run_errors_print_under_the_contract() {
        let plain = anyhow::anyhow!("inner").context("outer");
        assert_eq!(CliError::from_run(&plain), "ERROR: [run] outer: inner");
        let refused: anyhow::Error = CliError::new("output", "no").into();
        assert_eq!(CliError::from_run(&refused), "ERROR: [output] no");
    }

    /// The plugins directory is `--dir_plugins`, else `<dir>/Plugins`
    /// (Runner.pm:907), a leading `~` expanded as `glob` expands it (`:908`).
    #[test]
    fn dir_plugins_defaults_to_dir_slash_plugins() {
        let home = std::env::var("HOME").unwrap();
        assert_eq!(config(&["vep"]).dir_plugins, format!("{home}/.vep/Plugins"));
        assert_eq!(
            config(&["vep", "--dir", "/data/vep"]).dir_plugins,
            "/data/vep/Plugins"
        );
        assert_eq!(
            config(&["vep", "--dir", "/data/vep", "--dir_plugins", "/p"]).dir_plugins,
            "/p"
        );
        assert_eq!(
            config(&["vep", "--dir_plugins", "~/plugins"]).dir_plugins,
            format!("{home}/plugins")
        );
        assert!(!config(&["vep"]).safe);
        assert!(config(&["vep", "--safe"]).safe);
    }

    /// The warnings file is `--warning_file`, else `<output_file>_warnings.txt`
    /// with the output file as the user spelled it, `STDOUT` included
    /// (BaseVEP.pm:890 reads the `output_file` parameter as given).
    #[test]
    fn warning_file_defaults_beside_the_output_as_spelled() {
        assert_eq!(
            config(&["vep"]).warning_file,
            "variant_effect_output.txt_warnings.txt"
        );
        assert_eq!(
            config(&["vep", "-o", "out/x.txt"]).warning_file,
            "out/x.txt_warnings.txt"
        );
        assert_eq!(
            config(&["vep", "-o", "stdout"]).warning_file,
            "stdout_warnings.txt"
        );
        assert_eq!(
            config(&["vep", "-o", "a", "-o", "b"]).warning_file,
            "b_warnings.txt"
        );
        assert_eq!(
            config(&["vep", "-o", "x", "--warning_file", "w.txt"]).warning_file,
            "w.txt"
        );
        assert_eq!(
            config(&["vep", "--skipped_variants_file", "s.txt"])
                .skipped_variants_file
                .as_deref(),
            Some("s.txt")
        );
    }
}
