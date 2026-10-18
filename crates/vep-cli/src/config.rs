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

/// Environment variables vep-rs reads that Ensembl VEP has no counterpart for.
pub const ENVIRONMENT: &[(&str, &str)] = &[(
    "VEP_FASTA_STORAGE",
    "how the FASTA is held: `mmap` (default) maps the file, `memory` reads it into RAM; any other value is an error",
)];

/// The components an `ERROR: [<component>]` line names, with their meaning.
pub const COMPONENTS: &[(&str, &str)] = &[
    (
        "cli",
        "an unsupported flag, or a combination the option rules refuse",
    ),
    ("output", "the output file or format"),
    ("cache", "the transcript cache: its layout or its metadata"),
    (
        "filter",
        "a `--transcript_filter` Ensembl VEP would refuse, or whose regular expression is outside the supported dialect",
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
        state.apply_option_sets();
        state
            .check(state.is_active("safe"))
            .map_err(CliError::cli)?;
        for rule in OWN_RULES {
            if (rule.applies)(&state) {
                return Err(CliError::new(rule.component, rule.message));
            }
        }
        Ok(Config::resolve(args, &state))
    }

    /// The configuration from the arguments and the flag state after the option
    /// sets: a value an option set assigns comes from the state, everything else
    /// from the arguments.
    fn resolve(args: Args, state: &FlagState) -> Self {
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
        let output_file = if last.eq_ignore_ascii_case("STDOUT") {
            "STDOUT".to_string()
        } else {
            last
        };

        let dir = expand_tilde(&args.dir);
        let dir_cache = args
            .dir_cache
            .map(|d| expand_tilde(&d))
            .unwrap_or_else(|| dir.clone());

        let cache_version = args.cache_version.unwrap_or(vep_core::VEP_VERSION);

        let distance = parse_distance(&args.distance);

        let fields = args.fields.map(|f| {
            f.split(',')
                .map(|s| s.trim().to_string())
                .collect::<Vec<_>>()
        });

        if let Some(notice) = regulatory_unimplemented_notice(args.regulatory, args.everything) {
            tracing::warn!("{}", notice);
        }

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
            cache: on("cache"),
            offline: on("offline"),
            dir,
            dir_cache,
            dir_plugins: args.dir_plugins,
            cache_version,
            fasta: args.fasta,
            shift_3prime: args.shift_3prime,
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
            sift: state.scalar("sift").map(str::to_string),
            polyphen: state.scalar("polyphen").map(str::to_string),
            numbers: on("numbers") || vcf_output,
            domains: on("domains"),
            regulatory: on("regulatory"),
            variant_class: on("variant_class"),
            af: on("af"),
            af_1kg: on("af_1kg"),
            af_gnomade: on("af_gnomade"),
            af_gnomadg: on("af_gnomadg"),
            max_af: on("max_af"),
            pubmed: on("pubmed"),
            gene_phenotype: on("gene_phenotype"),
            mirna: on("mirna"),
            check_existing: on("check_existing"),
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
            plugin: args.plugin,
            custom: state.list("custom"),
            refseq: args.refseq,
            merged: args.merged,
            gencode_basic: args.gencode_basic,
            vcf_info_field: args.vcf_info_field,
            json_cache: args.json_cache,
            fields,
            dont_skip: args.dont_skip,
            allow_non_variant: on("allow_non_variant"),
        }
    }

    /// The cache directory the output headers name: the JSON cache when one is
    /// given, else Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout
    /// (CacheDir.pm:315 and :368), the species directory suffixed `_refseq` or
    /// `_merged` under those flags (:312), `<version>` alone when no assembly is set,
    /// consecutive `/` collapsed to one (:373), so a `--dir_cache` ending in `/` adds none.
    pub fn cache_dir_path(&self) -> String {
        if let Some(json_cache) = &self.json_cache {
            return json_cache.clone();
        }
        let mut species_dir = self.species.clone();
        if self.refseq {
            species_dir.push_str("_refseq");
        }
        if self.merged {
            species_dir.push_str("_merged");
        }
        let assembly_suffix = self
            .assembly
            .as_deref()
            .map(|a| format!("_{}", a))
            .unwrap_or_default();
        let joined = format!(
            "{}/{}/{}{}",
            self.dir_cache, species_dir, self.cache_version, assembly_suffix
        );
        let mut path = String::with_capacity(joined.len());
        for c in joined.chars() {
            if c == '/' && path.ends_with('/') {
                continue;
            }
            path.push(c);
        }
        path
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

    /// `--minimal` switches on `--allele_number` (the option set at Config.pm:524-529).
    #[test]
    fn minimal_implies_allele_number() {
        assert!(config(&["vep", "--minimal"]).allele_number);
        assert!(!config(&["vep"]).allele_number);
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

    /// Ensembl VEP's `CacheDir.pm` appends the species directory (:315), suffixed
    /// `_refseq` or `_merged` under those flags (:312), and then
    /// `<version>_<assembly>` (:368), and collapses consecutive `/` (:373), so
    /// `--dir_cache /opt/vep --assembly GRCh38` names
    /// `/opt/vep/homo_sapiens/<version>_GRCh38`, as do `--dir_cache /opt/vep/` and
    /// `--dir_cache /opt//vep`; `--refseq` names
    /// `/opt/vep/homo_sapiens_refseq/<version>_GRCh38` and `--merged`
    /// `/opt/vep/homo_sapiens_merged/<version>_GRCh38`.
    #[test]
    fn test_cache_dir_path() {
        let path = |flags: &[&str]| {
            let mut argv = vec!["vep", "--dir_cache", "/opt/vep", "--assembly", "GRCh38"];
            argv.extend_from_slice(flags);
            config(&argv).cache_dir_path()
        };
        let version = vep_core::VEP_VERSION;
        assert_eq!(path(&[]), format!("/opt/vep/homo_sapiens/{version}_GRCh38"));
        assert_eq!(
            path(&["--refseq"]),
            format!("/opt/vep/homo_sapiens_refseq/{version}_GRCh38")
        );
        assert_eq!(
            path(&["--merged"]),
            format!("/opt/vep/homo_sapiens_merged/{version}_GRCh38")
        );
        for dir_cache in ["/opt/vep/", "/opt//vep"] {
            assert_eq!(
                config(&["vep", "--dir_cache", dir_cache, "--assembly", "GRCh38"]).cache_dir_path(),
                format!("/opt/vep/homo_sapiens/{version}_GRCh38"),
                "{dir_cache}"
            );
        }
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
}
