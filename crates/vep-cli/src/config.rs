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
use std::str::FromStr;

use clap::FromArgMatches;

use vep_core::terms::Vocabulary;
use vep_io::input::Delimiter;

use crate::args::Args;
use crate::input_shape::InputShape;
use crate::nearest::NearestType;
use crate::option_rules::{split_list_flag, FlagState, OWN_RULES};

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
        "--gp is not supported: it reads a GRCh37 chr:pos from the INFO/GP field that some 1000 Genomes VCFs carry",
        "Config.pm:120; Parser/VCF.pm:271-285 takes chr:pos from INFO/GP and skips lines without it",
        "the record's own coordinates",
    ),
    refused(
        "ambiguous_hgvs",
        "--ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g., m.), transcript (c., n.) and protein (p.) notation, each resolved to one variant, a protein change with several equally short nucleotide paths skipped with Ensembl VEP's message",
        "Config.pm:202; Parser/HGVS.pm:178-186, 224 resolves the protein notation to every variant through the database; VariationFeatureAdaptor.pm:1956 refuses an ambiguous one without it",
        "genomic, transcript and protein HGVS input, one variant per notation",
    ),
    ErrorTableEntry {
        flag: "af_exac",
        declared: false,
        treatment: Treatment::Refused {
            component: "cli",
            message: "--af_exac was removed in Ensembl VEP release 107: use --af_gnomade",
        },
        ensembl_vep: "Config.pm:618 names it among --database's incompatibilities while @VEP_PARAMS (:73-276) does not declare it, so Ensembl VEP refuses it as an unexpected parameter (vep:30)",
        route: "--af_gnomade",
    },
    refused(
        "regulatory_gff",
        "--regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache; regulatory features are read from the cache alone",
        "Config.pm:210; the option set at :434-439 switches on --regulatory",
        "none",
    ),
    refused(
        "extended_promoters",
        "--extended_promoters is not supported: it would widen the promoters read with --regulatory_gff to their extended bounds; regulatory features are read from the cache alone",
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
        "the default it names, the input's reference allele, stands without it, and given beside --use_transcript_ref it keeps that default",
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
        "input",
        "the input: a file that does not exist, standard input without --format, a first data line of no known format, BCF input, CAid input (looked up in the Ensembl variation database, which vep-rs does not reach), an identifier list against a cache without variation files, a region, SPDI or HGVS input without --fasta, or a region input with --check_ref",
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
        "the cache: nothing resolves at the configured path, its layout or metadata is unusable, a transcript shard cannot be read (a transcript an HGVS input names included), its assembly disagrees with --assembly or with its own directory name, or a variation file cannot be opened or read",
    ),
    (
        "fasta",
        "the reference FASTA: the file or directory named cannot be read as an indexed FASTA, or a flag that reads it (`--hgvs`, `--hgvsg`, `--check_ref`, `--lookup_ref`, `--use_transcript_ref`, `--shift_3prime 1`) was given and none was named or found in the cache directory",
    ),
    (
        "filter",
        "a `--transcript_filter` Ensembl VEP would refuse, or whose regular expression is outside the supported dialect",
    ),
    (
        "check_svs",
        "the `--check_svs` file: it is not a file, has no `.csi` or `.tbi` index beside it, or cannot be read",
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

/// The `--input_data` string the run reads in place of its input file: the
/// flag's value when it is true in Perl (Runner.pm:818 reads it under `if`), so
/// an empty string or `0` is no input data and the input file stands.
pub fn input_data_for(args: &Args) -> Option<&str> {
    args.input_data
        .as_deref()
        .filter(|data| !data.is_empty() && *data != "0")
}

/// Where a command line's warnings go: `--warning_file` when its value is true
/// in Perl, else `<output_file>_warnings.txt` beside the last `-o` value as the
/// user spelled it, or beside `vep` when that value is false in Perl
/// (BaseVEP.pm:890 `$file ||= ($self->param('output_file') || 'vep')`, so an
/// empty value or `0` names no file on either side); readable before the
/// configuration is resolved, so the sink can take the warnings the resolution
/// itself raises.
pub fn warning_file_for(args: &Args) -> String {
    match args
        .warning_file
        .as_deref()
        .filter(|file| crate::argv::truthy(file))
    {
        Some(file) => file.to_string(),
        None => {
            let output_file = match args.output_file.last().map(String::as_str) {
                None => "variant_effect_output.txt",
                Some(file) if crate::argv::truthy(file) => file,
                Some(_) => "vep",
            };
            crate::warnings::default_warning_file(output_file)
        }
    }
}

/// A parsed command line: the argument struct and the flag state clap saw.
pub struct ParsedArgs {
    pub args: Args,
    pub state: FlagState,
    /// The tokens after the program name, as given, plus one `--<flag>=<value>`
    /// per switch a source set to a false value ([`Config::command_line`]).
    pub command_line: Vec<String>,
}

impl ParsedArgs {
    /// Records the switches a source behind the command line set to a false value
    /// ([`crate::argv::Prepassed::defined_off`], each with that value) as present
    /// but off, as Ensembl VEP's configuration holds them, and in the recorded
    /// command line as `--<flag>=<value>`, so the header prints `--<flag> <value>`
    /// as `full_command` prints a stored value other than `1` (Config.pm:745); a
    /// flag the command line gave keeps its value.
    pub fn with_defined_off<'a>(
        mut self,
        switches: impl IntoIterator<Item = (&'a str, &'a str)>,
    ) -> Self {
        for (flag, value) in switches {
            if !self.state.is_given(flag) {
                self.state.set(flag, value);
                self.command_line.push(format!("--{flag}={value}"));
            }
        }
        self
    }
}

/// Parses a command line (the binary name first). A bare `--help` or
/// `--version` wins over everything else on the line; then a flag the error
/// table refuses before parsing is reported; then clap's own syntax errors; the
/// result carries the arguments and the flag state the option rules act on.
pub fn parse_command_line<I, T>(argv: I) -> Result<ParsedArgs, ParseFailure>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    let mut argv = argv.into_iter().map(Into::into);
    let program = argv.next();
    let mut options: Vec<OsString> = argv.collect();
    let command_line: Vec<String> = options
        .iter()
        .map(|token| token.to_string_lossy().into_owned())
        .collect();
    // Getopt::Long runs with `pass_through` (Config.pm:62), which leaves a bare
    // `--` and every token after it in `@ARGV`, and vep refuses whatever is left
    // (vep:30) unless `--help` was read ahead of the `--` (:29); `--version`, this
    // engine's own flag, wins the same way.
    let separator = options
        .iter()
        .position(|a| a == "--")
        .unwrap_or(options.len());
    // The usage wins over every other token (vep:29 before vep:30): `pass_through`
    // leaves an unknown flag, a bare word, a refused value and a value-less option
    // in `@ARGV`, and `help` is tested first, so clap is handed the one token,
    // the help ahead of the version. `--help=1` is no such token: it is refused,
    // as vep:30 refuses it.
    let usage_token = ["--help", "-h", "--version", "-V"]
        .into_iter()
        .find(|usage| options[..separator].iter().any(|a| a == usage));
    if let Some(usage) = usage_token {
        options = vec![OsString::from(usage)];
    } else if separator < options.len() {
        return Err(ParseFailure::Refused(CliError::cli(format!(
            "Unexpected extra command-line parameter(s): {}",
            command_line[separator..].join(", ")
        ))));
    }
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
        if options.iter().any(names_it) {
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
    Ok(ParsedArgs {
        args,
        state,
        command_line,
    })
}

/// Fully resolved configuration for a VEP run.
#[derive(Debug)]
#[allow(dead_code)] // Many fields are parsed from CLI args but only consumed by specific code paths
pub struct Config {
    pub input_file: String,
    pub input_data: Option<String>,
    pub format: String,
    /// The column delimiter of an Ensembl-format or VCF input under `--format
    /// ensembl` or `--format vcf`; under `guess` the detected one is used instead
    /// (Parser.pm:387-415), and for every other format this holds the tab default.
    pub delimiter: Delimiter,

    /// The last `-o` value, `STDOUT` for standard output in any spelling.
    pub output_file: String,
    /// The last `-o` value as the user spelled it: the name BaseRunner.pm:203
    /// tests for an existing file before it compares the value with `STDOUT`
    /// (:213, :222), so a file named `STDOUT` in the working directory refuses the
    /// run as any other output file would.
    pub output_file_as_spelled: String,
    pub output_format: String,
    /// Only consulted when `output_format == "parquet"`.
    /// `"nested"` = one row per variant, CSQ subfields as `LIST<T>`.
    /// `"flat"` = one row per (variant × consequence) tuple.
    pub parquet_shape: String,
    /// Keep the Parquet intermediate TSV after DuckDB finalize (debug aid).
    pub keep_intermediate: bool,
    pub parquet_row_group_size: u64,
    pub compress_output: Option<String>,
    pub force_overwrite: bool,
    pub no_headers: bool,
    /// The vocabulary consequence terms print in (`--terms`).
    pub terms: Vocabulary,
    /// Leave HGVSp's `=` and a VCF gene symbol's whitespace unescaped: `--no_escape`,
    /// which `--json` implies through Ensembl VEP's `rest` option set; `--output_format
    /// json` implies it here as well (`docs/intended-divergences.md`).
    pub no_escape: bool,
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
    /// else a directory given to `--cache`, else `--dir`; an empty or `0`
    /// `--dir_cache` counts as absent (AnnotationSourceAdaptor.pm:127's `||`).
    pub dir_cache: String,
    /// The plugins directory: `--dir_plugins`, else `<dir>/Plugins`
    /// (Runner.pm:907), an empty or `0` value counting as absent, and a leading
    /// `~` expanded; any other pattern `glob` would expand (`:908`) is one
    /// literal path.
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
    /// 116 drops a record whose span, end minus start, is wider before annotation
    /// and writes nothing for it in any format (`Parser.pm` `validate_vf`), a
    /// sequence record included; vep-rs annotates it in the default, tab and
    /// Parquet outputs, carries its VCF line without consequences and omits it
    /// from the JSON output (`docs/intended-divergences.md`).
    pub max_sv_size: Option<u64>,
    pub fork: usize,
    /// BGZF decompression worker threads for a BGZF input in any format, as the
    /// flag was given (a plain gzip input decompresses on one thread);
    /// [`resolve_decompression_threads`] turns it into a count once the input's
    /// shape is known.
    pub decompression_threads: DecompressionThreads,
    pub distance: (u64, u64),
    /// How the batches move through the pipeline (`--pipeline_mode`).
    pub pipeline_mode: PipelineMode,
    /// Batches a pipeline stage may hold ahead of the next one
    /// (`--pipeline_depth`, resolved).
    pub pipeline_depth: usize,
    /// Input records per pool task; `None` is the built-in sizing
    /// (`--render_task_records 0`).
    pub render_task_records: Option<usize>,
    /// How rendered output reaches the file or standard output
    /// (`--writer_strategy`).
    pub writer_strategy: WriterStrategy,
    /// The `buffered` writer strategy's buffer, in bytes (`--output_buffer_size`,
    /// resolved).
    pub output_buffer_size: usize,
    /// Print the resolved runtime knobs and the input's shape, then exit.
    pub resolved_config: bool,
    /// The tokens the parser read after the program name, as given: the command
    /// line's own, and the flags of the `--config` file, `vep.ini` and the
    /// `VEP_*` variables that [`crate::argv::prepass`] laid out around them, plus
    /// a `--<flag>=<value>` for each switch a source set to a false value
    /// ([`ParsedArgs::with_defined_off`]). A token that is not UTF-8 is held with
    /// each invalid sequence replaced. The `## VEP command-line:` header line is
    /// rendered from this list.
    pub command_line: Vec<String>,

    /// `--everything` was given: the HGVS flags it switched on are dropped, not
    /// refused, when no FASTA is found (Runner.pm:745-754).
    pub everything: bool,
    pub symbol: bool,
    pub biotype: bool,
    pub canonical: bool,
    pub mane: bool,
    pub mane_select: bool,
    /// Set by `--mane`'s option set (Config.pm:580-586); no flag of its own.
    pub mane_plus_clinical: bool,
    pub tsl: bool,
    pub appris: bool,
    pub ccds: bool,
    pub protein: bool,
    /// Append the translation version to the ENSP identifier (`--protein_version`).
    pub protein_version: bool,
    pub uniprot: bool,
    pub xref_refseq: bool,
    /// Append the transcript version to the Feature column (`--transcript_version`).
    pub transcript_version: bool,
    /// Append the gene version to the Gene column (`--gene_version`).
    pub gene_version: bool,
    pub hgvs: bool,
    pub hgvsg: bool,
    /// `--hgvsg_use_accession`: HGVSg names the chromosome by the `NC_` synonym
    /// of the input name, when the synonyms table has one.
    pub hgvsg_use_accession: bool,
    /// `--spdi`: the SPDI column.
    pub spdi: bool,
    /// `--ga4gh_vrs`: the GA4GH_VRS column, refused with `--vcf`.
    pub ga4gh_vrs: bool,
    /// `--hgvsp_use_prediction`: HGVSp's description after `p.` is parenthesised.
    pub hgvsp_use_prediction: bool,
    /// `--remove_hgvsp_version`: HGVSp's reference name carries no version.
    pub remove_hgvsp_version: bool,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    /// `--humdiv`: the PolyPhen column reads the HumDiv matrix. Ensembl VEP keeps
    /// the choice as `polyphen_analysis`, `humvar` by default (Config.pm:310) and
    /// `humdiv` through the flag's option set (:517-522).
    pub humdiv: bool,
    pub numbers: bool,
    pub domains: bool,
    pub regulatory: bool,
    pub variant_class: bool,
    pub overlaps: bool,
    /// `--nearest <type>`: the NEAREST column names the nearest protein-coding
    /// transcription start sites by this identifier; `None` without the flag.
    pub nearest: Option<NearestType>,
    pub af: bool,
    pub af_1kg: bool,
    pub af_gnomade: bool,
    /// `--af_gnomad`: the older spelling of `--af_gnomade`, its own flag with
    /// the same fields (Config.pm:132); refused beside `--af_gnomade` (:619).
    pub af_gnomad: bool,
    pub af_gnomadg: bool,
    pub max_af: bool,
    pub pubmed: bool,
    pub old_maf: bool,
    /// `--check_frequency`, which `--filter_common` switches on (Config.pm:442-446):
    /// the frequency filter of `crate::frequency_filter` runs.
    pub check_frequency: bool,
    /// `--freq_pop`, as given (Ensembl VEP's default `1KG_ALL`, Config.pm:321).
    pub freq_pop: String,
    /// `--freq_freq`, the filter's threshold (default 0.01, Config.pm:319).
    pub freq_freq: f64,
    /// `--freq_gt_lt`: `gt` or `lt` (default `gt`, Config.pm:322).
    pub freq_gt_lt: String,
    /// `--freq_filter`: `exclude` or `include` (default `exclude`, Config.pm:320).
    pub freq_filter: String,
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
    /// `--check_svs <file>`: the indexed structural-variant VCF the `SV` column
    /// is read from; `None` without the flag.
    pub check_svs: Option<String>,
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
    /// The `--pick_order` criteria as given, `None` for the default order and an
    /// empty list for an empty value (`Config.pm:761` splits it to nothing).
    pub pick_order: Option<Vec<String>>,
    /// The `--transcript_filter` strings in command-line order, then the one
    /// `--exclude_predicted` appends (`Config.pm:566-571`); empty when neither
    /// flag was passed.
    pub transcript_filter: Vec<String>,
    pub allele_number: bool,
    pub show_ref_allele: bool,
    pub uploaded_allele: bool,
    pub ambiguity: bool,
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
    /// `--individual` as given, split on commas (`Config.pm:328`): `all` or the
    /// sample names; `None` without the flag. Switches `allow_non_variant` on
    /// (`:413-417`) and refuses `--minimal` and `--individual_zyg` (`:626`).
    pub individual: Option<Vec<String>>,
    /// `--individual_zyg` as given, split on commas (`Config.pm:329`): `all` or
    /// the sample names; `None` without the flag. Switches `allow_non_variant` on
    /// (`:420-424`).
    pub individual_zyg: Option<Vec<String>>,
    /// `--process_ref_homs` (`Config.pm:143`): every sequence record gains the
    /// reference allele's row, and a reference-homozygous individual annotates it.
    pub process_ref_homs: bool,
    /// `--chr` as given, split on commas (`Config.pm:759-761`); `None` without the
    /// flag. The names it admits are read by `validate::ChromosomeList`.
    pub chr: Option<Vec<String>>,
    /// `--check_ref`: the reference allele is checked against the FASTA.
    pub check_ref: bool,
    /// `--lookup_ref`: the reference allele is read from the FASTA.
    pub lookup_ref: bool,
    /// `--use_transcript_ref` (`Config.pm:254`): each row's reference allele is
    /// the transcript's, and `GIVEN_REF` and `USED_REF` are written; off again
    /// under `--use_given_ref` (`:588-590`).
    pub use_transcript_ref: bool,
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
/// The runner tests the flag with `defined`, so a `0` a configuration file or the
/// environment set it to ([`ParsedArgs::with_defined_off`]) switches
/// `check_existing` on too, while the synonyms themselves stay off.
fn var_synonyms_implies_check_existing(state: &mut FlagState) -> Option<&'static str> {
    if !state.is_given("var_synonyms") || state.is_active("check_existing") {
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
        let ParsedArgs {
            args,
            mut state,
            command_line,
        } = parsed;
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
        // Only the Ensembl and VCF parsers split on the delimiter (Parser/VEP_input.pm:90,
        // Parser/VCF.pm:129), and under `guess` the one detection chooses replaces the
        // flag's (Parser.pm:414), so the flag is read, and a bad value refused, only
        // when one of those two formats is named.
        let delimiter_is_read = matches!(state.scalar("format"), Some("ensembl" | "vcf"));
        let delimiter = match args.delimiter.as_deref() {
            Some(value) if delimiter_is_read => crate::delimiter::from_flag(value)
                .map_err(|message| CliError::new("config", message))?,
            _ => Delimiter::Tab,
        };
        // After every refusal: Ensembl VEP's `Runner.pm:117-120` runs once
        // `Config.pm:769` has checked the command line, so a refused one prints
        // its error alone.
        let notices = var_synonyms_implies_check_existing(&mut state)
            .into_iter()
            .collect();
        let mut config = Config::resolve(args, &state, distance, delimiter);
        config.notices = notices;
        config.command_line = command_line;
        Ok(config)
    }

    /// The configuration from the arguments and the flag state after the option
    /// sets: a value an option set assigns comes from the state, everything else
    /// from the arguments.
    fn resolve(args: Args, state: &FlagState, distance: (u64, u64), delimiter: Delimiter) -> Self {
        let on = |flag: &str| state.is_active(flag);
        let input_data = input_data_for(&args).map(str::to_string);
        let output_format = state.scalar("output_format").unwrap_or("vep").to_string();
        // The VCF CSQ layout names SYMBOL, BIOTYPE, EXON and INTRON whatever the
        // flags, so a VCF output fills them by either route; Ensembl VEP's `vcf`
        // option set does so for `--vcf` alone.
        let vcf_output = output_format == "vcf";
        let json_output = output_format == "json";

        // Getopt::Long accumulates `-o` (`output_file|o=s@`) and the vep script keeps
        // the last (vep:33); `STDOUT` is matched in any case (BaseRunner.pm:213).
        let output_file_as_spelled = args
            .output_file
            .last()
            .cloned()
            .unwrap_or_else(|| "variant_effect_output.txt".to_string());
        let warning_file = warning_file_for(&args);
        let output_file = if output_file_as_spelled.eq_ignore_ascii_case("STDOUT") {
            "STDOUT".to_string()
        } else {
            output_file_as_spelled.clone()
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
        // `dir_cache || dir` (AnnotationSourceAdaptor.pm:127) and
        // `dir_plugins || dir.'/Plugins'` (Runner.pm:907) fall through an empty
        // string and `0`, Perl's false values, to their defaults.
        let dir_cache = args
            .dir_cache
            .filter(|d| crate::argv::truthy(d))
            .map(|d| expand_tilde(&d))
            .or(cache_dir_value)
            .unwrap_or_else(|| dir.clone());
        let dir_plugins = args
            .dir_plugins
            .filter(|d| crate::argv::truthy(d))
            .map(|d| expand_tilde(&d))
            .unwrap_or_else(|| format!("{dir}/Plugins"));

        // `cache_version || db_version || software_version` (CacheDir.pm:323):
        // a `0` is Perl-false and falls through to the release.
        let cache_version = args
            .cache_version
            .filter(|version| *version != 0)
            .unwrap_or(vep_core::VEP_VERSION);

        let fields = args.fields.as_deref().map(split_list_flag);

        // `uc($file) eq 'STDIN'` (Parser.pm:235): the literal in any case.
        let input_file = if args.input_file.eq_ignore_ascii_case("STDIN") {
            "STDIN".to_string()
        } else {
            args.input_file
        };

        // `shift_hgvs` is a setting only as 0 or 1 (BaseRunner.pm:498), and 0 sets
        // `shift_3prime` to 0 whatever the command line passes (Runner.pm:792-795).
        let shift_hgvs = args.shift_hgvs.as_deref() != Some("0");
        let shift_3prime = shift_hgvs && args.shift_3prime.as_deref() == Some("1");

        Config {
            input_file,
            input_data,
            format: state.scalar("format").unwrap_or("guess").to_string(),
            delimiter,
            output_file,
            output_file_as_spelled,
            output_format,
            parquet_shape: args.parquet_shape,
            keep_intermediate: args.keep_intermediate,
            parquet_row_group_size: args.parquet_row_group_size,
            compress_output: args.compress_output,
            force_overwrite: args.force_overwrite,
            no_headers: args.no_headers,
            // clap's enumeration admits the `%VALID` spellings alone; a `None` here
            // means it drifted from `Vocabulary::from_flag`, and the default is safer
            // than a crash.
            terms: Vocabulary::from_flag(&args.terms).unwrap_or_default(),
            no_escape: on("no_escape") || json_output,
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
            // clap validated the text with `decompression_threads_text`.
            decompression_threads: args.decompression_threads.parse().unwrap_or_default(),
            distance,
            pipeline_mode: args.pipeline_mode.parse().unwrap_or_default(),
            pipeline_depth: resolve_pipeline_depth(args.pipeline_depth),
            render_task_records: resolve_render_task_records(args.render_task_records),
            writer_strategy: args.writer_strategy.parse().unwrap_or_default(),
            output_buffer_size: resolve_output_buffer_size(args.output_buffer_size),
            resolved_config: args.resolved_config,
            command_line: Vec::new(),
            everything: on("everything"),
            symbol: on("symbol") || vcf_output,
            biotype: on("biotype") || vcf_output,
            canonical: on("canonical"),
            mane: on("mane"),
            mane_select: on("mane_select"),
            mane_plus_clinical: on("mane_plus_clinical"),
            tsl: on("tsl"),
            appris: on("appris"),
            ccds: on("ccds"),
            protein: on("protein"),
            protein_version: args.protein_version,
            uniprot: on("uniprot"),
            xref_refseq: args.xref_refseq,
            transcript_version: args.transcript_version,
            gene_version: args.gene_version,
            hgvs: on("hgvs"),
            hgvsg: args.hgvsg,
            hgvsg_use_accession: args.hgvsg_use_accession,
            spdi: args.spdi,
            ga4gh_vrs: args.ga4gh_vrs,
            hgvsp_use_prediction: args.hgvsp_use_prediction,
            remove_hgvsp_version: args.remove_hgvsp_version,
            sift: state.scalar("sift").map(str::to_string),
            polyphen: state.scalar("polyphen").map(str::to_string),
            humdiv: state.scalar("polyphen_analysis") == Some("humdiv"),
            numbers: on("numbers") || vcf_output,
            domains: on("domains"),
            regulatory: on("regulatory"),
            variant_class: on("variant_class"),
            overlaps: args.overlaps,
            // A value here passed clap's list and `VALID`, so it parses.
            nearest: state.scalar("nearest").and_then(|v| v.parse().ok()),
            af: on("af"),
            af_1kg: on("af_1kg"),
            af_gnomade: on("af_gnomade"),
            af_gnomad: on("af_gnomad"),
            af_gnomadg: on("af_gnomadg"),
            max_af: on("max_af"),
            pubmed: on("pubmed"),
            old_maf: on("old_maf"),
            check_frequency: on("check_frequency"),
            freq_pop: args.freq_pop,
            freq_freq: args.freq_freq,
            freq_gt_lt: args.freq_gt_lt,
            freq_filter: args.freq_filter,
            gene_phenotype: on("gene_phenotype"),
            mirna: on("mirna"),
            check_existing: on("check_existing"),
            no_check_alleles: args.no_check_alleles,
            exclude_null_alleles: args.exclude_null_alleles,
            failed: args.failed,
            clin_sig_allele: args.clin_sig_allele != 0,
            clinvar_somatic_classification: args.clinvar_somatic_classification,
            var_synonyms: on("var_synonyms"),
            check_svs: args.check_svs.as_deref().map(expand_tilde),
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
            pick_order: args.pick_order.map(|order| order.0),
            transcript_filter: state.list("transcript_filter"),
            allele_number: on("allele_number"),
            show_ref_allele: args.show_ref_allele,
            uploaded_allele: args.uploaded_allele,
            ambiguity: args.ambiguity,
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
            individual: state
                .is_given("individual")
                .then(|| state.list("individual")),
            individual_zyg: state
                .is_given("individual_zyg")
                .then(|| state.list("individual_zyg")),
            process_ref_homs: args.process_ref_homs,
            chr: state.is_given("chr").then(|| state.list("chr")),
            check_ref: args.check_ref,
            lookup_ref: args.lookup_ref,
            use_transcript_ref: on("use_transcript_ref"),
        }
    }
}

/// Cap for the `0` rule of `--decompression_threads`, the default: decompression
/// throughput plateaus between 4 and 10 threads and is flat past that, so the
/// rule stops at 10 rather than oversubscribing large hosts.
pub const DECOMPRESSION_THREADS_AUTO_CAP: usize = 10;

/// Auto-detect cap for `--fork=0`. One reader thread parses and hands out the
/// batches, so the wall time stops falling past about 32 workers and rises on
/// wider hosts, where the extra workers wait at the batch barrier and contend
/// for the allocator. User-supplied values are not capped.
const FORK_AUTO_CAP: usize = 32;

/// Resolve the --fork value. `0` means "auto": use the number of available
/// logical CPUs, capped at `FORK_AUTO_CAP`. Any non-zero value is returned
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

/// An input whose records inflate to at least this many bytes is bound by
/// inflation rather than annotation; a wider inflater pool only pays there.
pub const INFLATE_BOUND_BYTES_PER_RECORD: u64 = 8_000;

/// Inflaters for an input bound by inflation on an aarch64 host, before the
/// host's parallelism caps them: the count past which such an input stopped
/// gaining.
#[cfg(target_arch = "aarch64")]
pub const INFLATE_BOUND_INFLATERS: usize = 32;

/// Inflaters for an input bound by inflation on any other architecture: on
/// x86_64 no count above the `0` rule's cap read such an input faster, so the
/// `auto` rule keeps that cap there.
#[cfg(not(target_arch = "aarch64"))]
pub const INFLATE_BOUND_INFLATERS: usize = DECOMPRESSION_THREADS_AUTO_CAP;

/// The `--decompression_threads` value: the form the flag was given in, resolved
/// to a count once the input's shape is known.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecompressionThreads {
    /// `0`: the host's available parallelism, capped at
    /// [`DECOMPRESSION_THREADS_AUTO_CAP`].
    #[default]
    HostCapped,
    /// `auto`: the count the input's shape asks for
    /// ([`resolve_decompression_threads_auto`]).
    Auto,
    /// A positive count, lowered to the host's available parallelism.
    Fixed(usize),
}

impl FromStr for DecompressionThreads {
    type Err = String;

    /// Accepts exactly `0`, `auto` and a positive integer written in digits
    /// alone (`usize` also parses a leading `+`, which the flag refuses).
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "0" {
            return Ok(DecompressionThreads::HostCapped);
        }
        if s == "auto" {
            return Ok(DecompressionThreads::Auto);
        }
        let digits_only = s.bytes().all(|b| b.is_ascii_digit());
        match s.parse::<usize>() {
            Ok(count) if count > 0 && digits_only => Ok(DecompressionThreads::Fixed(count)),
            _ => Err(format!(
                "'{s}' is not 0, auto or a positive whole number of threads"
            )),
        }
    }
}

impl fmt::Display for DecompressionThreads {
    /// The value as the flag spells it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DecompressionThreads::HostCapped => f.write_str("0"),
            DecompressionThreads::Auto => f.write_str("auto"),
            DecompressionThreads::Fixed(n) => write!(f, "{n}"),
        }
    }
}

/// The `--decompression_threads` text, validated as [`DecompressionThreads`]
/// parses it; the command-line parser's value parser.
pub fn decompression_threads_text(text: &str) -> Result<String, String> {
    text.parse::<DecompressionThreads>()
        .map(|_| text.to_string())
}

/// The logical CPUs the decompression rules assume on a host that reports no
/// count.
const UNREPORTED_HOST_CPUS: usize = 4;

/// The host's logical CPU count, `None` when the host does not report one.
pub fn reported_parallelism() -> Option<usize> {
    std::thread::available_parallelism().ok().map(|n| n.get())
}

/// The host's logical CPUs as the decompression rules read them:
/// `UNREPORTED_HOST_CPUS` when the host does not report a count.
pub fn host_parallelism() -> usize {
    reported_parallelism().unwrap_or(UNREPORTED_HOST_CPUS)
}

/// Resolve the --decompression_threads value against the input's shape, on a
/// host of [`host_parallelism`] CPUs: a fixed count lowered to the host's CPUs
/// when the host reports them (the block reader starts one thread and one block
/// buffer per worker, so a count far above them would exhaust threads or memory
/// before it read a block, and workers beyond the CPUs decompress nothing more),
/// else kept as given; `0` as the host's CPUs capped at
/// [`DECOMPRESSION_THREADS_AUTO_CAP`]; `auto` by
/// [`resolve_decompression_threads_auto`].
pub fn resolve_decompression_threads(value: DecompressionThreads, shape: &InputShape) -> usize {
    resolve_decompression_threads_on(value, shape, reported_parallelism())
}

/// [`resolve_decompression_threads`] on a host whose logical CPU count is
/// `reported`, `None` when the host does not report one.
fn resolve_decompression_threads_on(
    value: DecompressionThreads,
    shape: &InputShape,
    reported: Option<usize>,
) -> usize {
    let cpus = reported.unwrap_or(UNREPORTED_HOST_CPUS);
    match value {
        DecompressionThreads::Fixed(n) => reported.map_or(n, |reported_cpus| n.min(reported_cpus)),
        DecompressionThreads::HostCapped => cpus.min(DECOMPRESSION_THREADS_AUTO_CAP),
        DecompressionThreads::Auto => resolve_decompression_threads_auto(shape, cpus),
    }
}

/// The `auto` rule: an indexed input whose records inflate to at least
/// [`INFLATE_BOUND_BYTES_PER_RECORD`] bytes gets [`INFLATE_BOUND_INFLATERS`]
/// inflaters, capped at `available_parallelism`; any other shape gets the `0`
/// rule, `available_parallelism` capped at [`DECOMPRESSION_THREADS_AUTO_CAP`].
pub fn resolve_decompression_threads_auto(
    shape: &InputShape,
    available_parallelism: usize,
) -> usize {
    match shape {
        InputShape::Indexed {
            inflated_bytes_per_record,
            ..
        } if *inflated_bytes_per_record >= INFLATE_BOUND_BYTES_PER_RECORD => {
            INFLATE_BOUND_INFLATERS.min(available_parallelism)
        }
        _ => available_parallelism.min(DECOMPRESSION_THREADS_AUTO_CAP),
    }
}

/// How the batches move through the pipeline (`--pipeline_mode`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PipelineMode {
    /// Each batch is parsed, annotated and rendered on the worker pool as a
    /// whole before the next one starts.
    Batched,
    /// Each batch is handed to the pool as independent tasks that parse,
    /// annotate and render their own records; the next batch follows without
    /// waiting for them and the writer restores the input order.
    #[default]
    Streaming,
}

impl PipelineMode {
    /// The value as the flag spells it.
    pub fn name(self) -> &'static str {
        match self {
            PipelineMode::Batched => "batched",
            PipelineMode::Streaming => "streaming",
        }
    }
}

impl FromStr for PipelineMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "batched" => Ok(PipelineMode::Batched),
            "streaming" => Ok(PipelineMode::Streaming),
            other => Err(format!("unknown pipeline mode '{other}'")),
        }
    }
}

/// How rendered output reaches the file or standard output (`--writer_strategy`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WriterStrategy {
    /// Through a buffer of [`Config::output_buffer_size`] bytes.
    #[default]
    Buffered,
    /// Each rendered chunk with one write call.
    Direct,
    /// The chunks that are ready together with one vectored write call,
    /// repeated until every byte is written.
    Vectored,
}

impl WriterStrategy {
    /// The value as the flag spells it.
    pub fn name(self) -> &'static str {
        match self {
            WriterStrategy::Buffered => "buffered",
            WriterStrategy::Direct => "direct",
            WriterStrategy::Vectored => "vectored",
        }
    }
}

impl FromStr for WriterStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "buffered" => Ok(WriterStrategy::Buffered),
            "direct" => Ok(WriterStrategy::Direct),
            "vectored" => Ok(WriterStrategy::Vectored),
            other => Err(format!("unknown writer strategy '{other}'")),
        }
    }
}

/// Resolve the `--pipeline_depth` value: `0` is `crate::pipeline::PIPELINE_DEPTH`.
/// The flag's parser bounds the value, so it fits a `usize` on every target.
pub fn resolve_pipeline_depth(raw: u64) -> usize {
    if raw == 0 {
        crate::pipeline::PIPELINE_DEPTH
    } else {
        raw as usize
    }
}

/// Resolve the `--render_task_records` value: `0` is `None`, the built-in sizing
/// (per batch from the worker count on the fused path, a fixed chunk on the
/// plugin path).
pub fn resolve_render_task_records(raw: u64) -> Option<usize> {
    if raw == 0 {
        None
    } else {
        Some(raw as usize)
    }
}

/// Resolve the `--output_buffer_size` value: `0` is
/// `crate::runner::OUTPUT_BUFFER_BYTES`.
pub fn resolve_output_buffer_size(raw: u64) -> usize {
    if raw == 0 {
        crate::runner::OUTPUT_BUFFER_BYTES
    } else {
        raw as usize
    }
}

/// The `--distance` value as (upstream, downstream) bp: one integer for both, or
/// the pair `up,down`. The value is a `@LIST_FLAGS` entry, split with Perl's
/// `split` (`Config.pm:761`, [`split_list_flag`]), so `5000,` is the one integer
/// `5000`. Ensembl VEP reads the pair with `$d ||= $u` (`BaseRunner.pm` 506),
/// which fires on the false string `0` alone, so the downstream part `0` is the
/// upstream value while `00` or `+0` is a downstream distance of 0. Ensembl VEP
/// numifies a non-numeric part to 0 under an `Argument "x" isn't numeric` warning
/// and drops the up- or downstream rows that part governs (both for a single
/// value), and reads no part after the second; here each is a refusal naming the
/// value, as is a part padded with whitespace, which Perl numifies but `||=`
/// reads as set.
fn parse_distance(s: &str) -> Result<(u64, u64), String> {
    let refused = || format!("--distance expects one or two integers (up[,down]); got {s:?}");
    let integer = |part: &str| part.parse::<u64>().map_err(|_| refused());
    let parts = split_list_flag(s);
    let (up, down) = match parts.as_slice() {
        [up] => (up.as_str(), up.as_str()),
        [up, down] if down == "0" => (up.as_str(), up.as_str()),
        [up, down] => (up.as_str(), down.as_str()),
        _ => return Err(refused()),
    };
    Ok((integer(up)?, integer(down)?))
}

/// Expand a leading `~` to the user's home directory
/// ([`crate::argv::expand_tilde_text`] with `HOME`, empty when it is unset).
fn expand_tilde(path: &str) -> String {
    crate::argv::expand_tilde_text(path, &std::env::var("HOME").unwrap_or_default())
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

    /// The kind of clap's answer to a command line it does not parse into a run.
    fn usage_kind(argv: &[&str]) -> clap::error::ErrorKind {
        match parse_command_line(argv) {
            Err(ParseFailure::Usage(e)) => e.kind(),
            Err(ParseFailure::Refused(e)) => panic!("refused, not a usage error: {e}"),
            Ok(_) => panic!("{argv:?} parsed into a run"),
        }
    }

    /// A bare `--help` or `--version` wins over every other token of the line
    /// (vep:29 before vep:30 under `pass_through`, Config.pm:62): a flag the
    /// error table refuses, an unknown flag, a bare word, a refused value and an
    /// option left without its value are all dropped, the help ahead of the
    /// version; `--help=1` is no such token and is clap's refusal, as vep:30
    /// refuses it; behind `--` either is an extra parameter.
    #[test]
    fn help_and_version_win_over_every_other_token() {
        use clap::error::ErrorKind;
        for argv in [
            &["vep", "--af_exac", "--help"][..],
            &["vep", "--help", "--af_exac"][..],
            &["vep", "--bogus", "--help"][..],
            &["vep", "--stats_html", "--help"][..],
            &["vep", "extra", "--help"][..],
            &["vep", "--fork=x", "--help"][..],
            &["vep", "--fork", "--help"][..],
            &["vep", "--help", "--input_file"][..],
            &["vep", "--version", "--help"][..],
            &["vep", "-h", "--bogus"][..],
            &["vep", "--help", "--", "foo"][..],
        ] {
            assert_eq!(usage_kind(argv), ErrorKind::DisplayHelp, "{argv:?}");
        }
        for argv in [
            &["vep", "--bogus", "--version"][..],
            &["vep", "--af_exac", "-V"][..],
            &["vep", "--fork", "--version"][..],
        ] {
            assert_eq!(usage_kind(argv), ErrorKind::DisplayVersion, "{argv:?}");
        }
        let kind = usage_kind(&["vep", "--help=1"]);
        assert!(
            kind != ErrorKind::DisplayHelp && kind != ErrorKind::DisplayVersion,
            "{kind:?}"
        );
        assert_eq!(
            refusal(&["vep", "--", "--help"]),
            "ERROR: [cli] Unexpected extra command-line parameter(s): --, --help"
        );
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

    /// `0` on each numeric runtime knob is the built-in value; any other value
    /// passes through.
    #[test]
    fn zero_resolves_to_the_built_in_runtime_values() {
        assert_eq!(resolve_pipeline_depth(0), crate::pipeline::PIPELINE_DEPTH);
        assert_eq!(resolve_pipeline_depth(1), 1);
        assert_eq!(resolve_pipeline_depth(8), 8);
        assert_eq!(resolve_render_task_records(0), None);
        assert_eq!(resolve_render_task_records(1), Some(1));
        assert_eq!(resolve_render_task_records(64), Some(64));
        assert_eq!(
            resolve_output_buffer_size(0),
            crate::runner::OUTPUT_BUFFER_BYTES
        );
        assert_eq!(resolve_output_buffer_size(1), 1);
        assert_eq!(resolve_output_buffer_size(4096), 4096);
    }

    /// `--decompression_threads` takes exactly `0`, `auto` and a positive whole
    /// number; everything else is clap's usage error.
    #[test]
    fn decompression_threads_parses_its_three_forms() {
        for (text, value) in [
            ("0", DecompressionThreads::HostCapped),
            ("auto", DecompressionThreads::Auto),
            ("1", DecompressionThreads::Fixed(1)),
            ("32", DecompressionThreads::Fixed(32)),
        ] {
            assert_eq!(
                text.parse::<DecompressionThreads>().unwrap(),
                value,
                "{text}"
            );
            assert_eq!(decompression_threads_text(text).unwrap(), text);
            assert_eq!(
                config(&["vep", "--decompression_threads", text]).decompression_threads,
                value,
                "{text}"
            );
        }
        assert_eq!(
            config(&["vep"]).decompression_threads,
            DecompressionThreads::HostCapped
        );
        assert_eq!(DecompressionThreads::Fixed(7).to_string(), "7");
        assert_eq!(DecompressionThreads::Auto.to_string(), "auto");
        assert_eq!(DecompressionThreads::HostCapped.to_string(), "0");
        for text in ["-1", "1.5", "AUTO", "", "+1", "0x10", "auto1"] {
            assert!(text.parse::<DecompressionThreads>().is_err(), "{text}");
            assert!(decompression_threads_text(text).is_err(), "{text}");
            match Config::from_argv(["vep", "--decompression_threads", text]).unwrap_err() {
                ParseFailure::Usage(e) => assert!(
                    matches!(
                        e.kind(),
                        clap::error::ErrorKind::ValueValidation
                            | clap::error::ErrorKind::UnknownArgument
                    ),
                    "{text}: {e}"
                ),
                ParseFailure::Refused(e) => panic!("{text} refused instead: {e}"),
            }
        }
    }

    /// The `auto` rule over the three shapes and both sides of the threshold:
    /// an unknown shape and an input below the threshold take the `0` rule, an
    /// input at or above it takes the inflate-bound pool, and the host's
    /// parallelism caps both.
    #[test]
    fn auto_decompression_threads_follow_the_input_shape() {
        let indexed = |inflated_bytes_per_record: u64| InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate: 1_000,
            file_bytes: 1 << 20,
        };
        let wide = 64;
        assert_eq!(
            resolve_decompression_threads_auto(&InputShape::Unknown, wide),
            DECOMPRESSION_THREADS_AUTO_CAP
        );
        assert_eq!(
            resolve_decompression_threads_auto(&InputShape::Unknown, 4),
            4
        );
        assert_eq!(
            resolve_decompression_threads_auto(&indexed(INFLATE_BOUND_BYTES_PER_RECORD - 1), wide),
            DECOMPRESSION_THREADS_AUTO_CAP
        );
        assert_eq!(
            resolve_decompression_threads_auto(&indexed(INFLATE_BOUND_BYTES_PER_RECORD), wide),
            INFLATE_BOUND_INFLATERS
        );
        assert_eq!(
            resolve_decompression_threads_auto(&indexed(INFLATE_BOUND_BYTES_PER_RECORD * 2), 16),
            INFLATE_BOUND_INFLATERS.min(16)
        );
        assert_eq!(
            resolve_decompression_threads_auto(&indexed(INFLATE_BOUND_BYTES_PER_RECORD * 2), 1),
            1
        );

        let host = host_parallelism();
        assert_eq!(
            resolve_decompression_threads(DecompressionThreads::Fixed(3), &InputShape::Unknown),
            3.min(host)
        );
        assert_eq!(
            resolve_decompression_threads(DecompressionThreads::HostCapped, &indexed(1 << 20)),
            host.min(DECOMPRESSION_THREADS_AUTO_CAP)
        );
        assert_eq!(
            resolve_decompression_threads(DecompressionThreads::Auto, &indexed(1 << 20)),
            INFLATE_BOUND_INFLATERS.min(host)
        );
        assert_eq!(
            resolve_decompression_threads(DecompressionThreads::Auto, &InputShape::Unknown),
            host.min(DECOMPRESSION_THREADS_AUTO_CAP)
        );
    }

    /// Every value of each runtime knob reaches its resolved `Config` field, an
    /// absent flag resolves as `0` does and the modes and strategies as their
    /// defaults (`streaming`, `buffered`), and a value outside the enumeration or
    /// the bounds is clap's usage error.
    #[test]
    fn runtime_knobs_reach_the_config() {
        let plain = config(&["vep"]);
        assert!(!plain.resolved_config);
        assert!(config(&["vep", "--resolved_config"]).resolved_config);
        assert_eq!(plain.pipeline_mode, PipelineMode::Streaming);
        assert_eq!(plain.pipeline_mode, PipelineMode::default());
        assert_eq!(plain.pipeline_depth, crate::pipeline::PIPELINE_DEPTH);
        assert_eq!(plain.render_task_records, None);
        assert_eq!(plain.writer_strategy, WriterStrategy::Buffered);
        assert_eq!(plain.writer_strategy, WriterStrategy::default());
        assert_eq!(plain.output_buffer_size, crate::runner::OUTPUT_BUFFER_BYTES);

        for (value, mode) in [
            ("batched", PipelineMode::Batched),
            ("streaming", PipelineMode::Streaming),
        ] {
            let c = config(&["vep", "--pipeline_mode", value]);
            assert_eq!(c.pipeline_mode, mode, "{value}");
            assert_eq!(mode.name(), value);
            assert_eq!(value.parse::<PipelineMode>().unwrap(), mode);
        }
        for (value, strategy) in [
            ("buffered", WriterStrategy::Buffered),
            ("direct", WriterStrategy::Direct),
            ("vectored", WriterStrategy::Vectored),
        ] {
            let c = config(&["vep", "--writer_strategy", value]);
            assert_eq!(c.writer_strategy, strategy, "{value}");
            assert_eq!(strategy.name(), value);
            assert_eq!(value.parse::<WriterStrategy>().unwrap(), strategy);
        }
        assert!("Batched".parse::<PipelineMode>().is_err());
        assert!("buffer".parse::<WriterStrategy>().is_err());

        let zeros = config(&[
            "vep",
            "--pipeline_depth",
            "0",
            "--render_task_records",
            "0",
            "--output_buffer_size",
            "0",
        ]);
        assert_eq!(zeros.pipeline_depth, crate::pipeline::PIPELINE_DEPTH);
        assert_eq!(zeros.render_task_records, None);
        assert_eq!(zeros.output_buffer_size, crate::runner::OUTPUT_BUFFER_BYTES);
        let set = config(&[
            "vep",
            "--pipeline_depth",
            "8",
            "--render_task_records",
            "7",
            "--output_buffer_size",
            "4096",
        ]);
        assert_eq!(set.pipeline_depth, 8);
        assert_eq!(set.render_task_records, Some(7));
        assert_eq!(set.output_buffer_size, 4096);
        let bounds = config(&[
            "vep",
            "--pipeline_depth",
            "1024",
            "--render_task_records",
            "1000000",
            "--output_buffer_size",
            "1073741824",
        ]);
        assert_eq!(bounds.pipeline_depth, 1024);
        assert_eq!(bounds.render_task_records, Some(1_000_000));
        assert_eq!(bounds.output_buffer_size, 1_073_741_824);

        for argv in [
            &["vep", "--pipeline_mode", "fused"][..],
            &["vep", "--writer_strategy", "mmap"][..],
            &["vep", "--pipeline_depth", "-1"][..],
            &["vep", "--pipeline_depth", "1025"][..],
            &["vep", "--pipeline_depth", "18446744073709551615"][..],
            &["vep", "--render_task_records", "many"][..],
            &["vep", "--render_task_records", "1000001"][..],
            &["vep", "--output_buffer_size", "4k"][..],
            &["vep", "--output_buffer_size", "1073741825"][..],
            &["vep", "--output_buffer_size", "18446744073709551615"][..],
        ] {
            match Config::from_argv(argv).unwrap_err() {
                ParseFailure::Usage(e) => assert!(
                    matches!(
                        e.kind(),
                        clap::error::ErrorKind::InvalidValue
                            | clap::error::ErrorKind::ValueValidation
                            | clap::error::ErrorKind::UnknownArgument
                    ),
                    "{argv:?}: {e}"
                ),
                ParseFailure::Refused(e) => panic!("{argv:?} refused instead: {e}"),
            }
        }
    }

    /// The `0` decompression worker count stays within its cap; a pinned count
    /// is kept up to the host's logical CPUs and lowered to them above that,
    /// and kept whole on a host that reports no count.
    #[test]
    fn decompression_threads_pinned_up_to_the_host_cpus() {
        let cpus = std::thread::available_parallelism().unwrap().get();
        let resolve = |value| resolve_decompression_threads(value, &InputShape::Unknown);
        let host_capped = resolve(DecompressionThreads::HostCapped);
        assert!((1..=DECOMPRESSION_THREADS_AUTO_CAP).contains(&host_capped));
        assert_eq!(resolve(DecompressionThreads::Fixed(1)), 1);
        assert_eq!(resolve(DecompressionThreads::Fixed(cpus)), cpus);
        assert_eq!(resolve(DecompressionThreads::Fixed(cpus + 1)), cpus);
        assert_eq!(resolve(DecompressionThreads::Fixed(usize::MAX)), cpus);

        let unreported =
            |value| resolve_decompression_threads_on(value, &InputShape::Unknown, None);
        assert_eq!(unreported(DecompressionThreads::Fixed(1)), 1);
        assert_eq!(unreported(DecompressionThreads::Fixed(64)), 64);
        assert_eq!(
            unreported(DecompressionThreads::HostCapped),
            UNREPORTED_HOST_CPUS
        );
        assert_eq!(unreported(DecompressionThreads::Auto), UNREPORTED_HOST_CPUS);
        assert_eq!(
            resolve_decompression_threads_on(
                DecompressionThreads::Fixed(64),
                &InputShape::Unknown,
                Some(16)
            ),
            16
        );
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
        let flags = [
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
            ("af_gnomad", c.af_gnomad),
            ("af_gnomadg", c.af_gnomadg),
            ("max_af", c.max_af),
            ("pubmed", c.pubmed),
            ("uniprot", c.uniprot),
            ("mane", c.mane),
            ("mane_select", c.mane_select),
            ("mane_plus_clinical", c.mane_plus_clinical),
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
            ("transcript_version", c.transcript_version),
            ("gene_version", c.gene_version),
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

    /// The `--everything` help, which the CLI reference prints, states what its
    /// `--regulatory` part writes, the regulatory and motif rows, and the two
    /// omissions, motif scores and cell-type activity.
    #[test]
    fn everything_help_names_the_regulatory_rows_and_omissions() {
        let mut command = Args::command_with_error_table();
        command.build();
        let help = command
            .get_arguments()
            .find(|a| a.get_long() == Some("everything"))
            .and_then(|a| a.get_long_help().or(a.get_help()))
            .map(|h| h.to_string())
            .unwrap();
        for phrase in [
            "regulatory and motif rows",
            "motif scores",
            "cell-type activity",
            "omitted",
        ] {
            assert!(help.contains(phrase), "{phrase:?} missing from {help}");
        }
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

    /// `--terms` takes Ensembl VEP's `%VALID` values (Config.pm:598: `SO`,
    /// `display`, `NCBI`), `SO` by default (:307), with `-t` as its short form
    /// (:176 `terms|t=s`); a value outside the list, `ensembl` included, is a
    /// value outside the enumeration and so clap's own error (exit 2), as for
    /// every declared `%VALID` flag.
    #[test]
    fn terms_rejects_ensembl_as_a_value_outside_the_enumeration() {
        assert_eq!(config(&["vep"]).terms, Vocabulary::So);
        assert_eq!(
            config(&["vep", "--terms", "display"]).terms,
            Vocabulary::Display
        );
        assert_eq!(config(&["vep", "-t", "NCBI"]).terms, Vocabulary::Ncbi);
        for bad in ["ensembl", "so", "ncbi"] {
            match Config::from_argv(["vep", "--terms", bad]).unwrap_err() {
                ParseFailure::Usage(e) => {
                    assert_eq!(e.kind(), clap::error::ErrorKind::InvalidValue, "{bad}")
                }
                ParseFailure::Refused(e) => panic!("{bad} refused instead: {e}"),
            }
        }
    }

    /// `--no_escape` is a flag of its own (Config.pm:228); `--json` implies it
    /// through the `rest` option set (:484-498), and `--output_format json` is
    /// read the same way here, where Ensembl VEP escapes that output.
    #[test]
    fn json_flag_implies_no_escape() {
        assert!(config(&["vep", "--no_escape"]).no_escape);
        assert!(config(&["vep", "--json"]).no_escape);
        assert!(config(&["vep", "--output_format", "json"]).no_escape);
        assert!(!config(&["vep"]).no_escape);
        assert!(!config(&["vep", "--vcf"]).no_escape);
        assert!(!config(&["vep", "--tab"]).no_escape);
    }

    /// `--uploaded_allele` (Config.pm:234) and `--ambiguity` (:229) are plain
    /// switches that no option set implies.
    #[test]
    fn uploaded_allele_and_ambiguity_are_plain_switches() {
        let c = config(&["vep", "--uploaded_allele", "--ambiguity"]);
        assert!(c.uploaded_allele && c.ambiguity);
        let plain = config(&["vep", "--everything", "--minimal", "--vcf"]);
        assert!(!plain.uploaded_allele && !plain.ambiguity);
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

        // A `0` or an empty value from a configuration file or the environment:
        // defined but false (Runner.pm:117 tests `defined`), so the implication
        // runs and the flag itself stays off, and the recorded command line
        // carries the value the header prints (Config.pm:745); a flag the
        // command line gave keeps its value and its line.
        for value in ["0", ""] {
            let parsed = parse_command_line(["vep"])
                .unwrap()
                .with_defined_off([("var_synonyms", value)]);
            assert_eq!(parsed.command_line, [format!("--var_synonyms={value}")]);
            let defined_off = Config::from_parsed(parsed).unwrap();
            assert!(!defined_off.var_synonyms && defined_off.check_existing);
            assert_eq!(defined_off.notices, [NOTICE]);
        }
        let given = parse_command_line(["vep", "--var_synonyms"])
            .unwrap()
            .with_defined_off([("var_synonyms", "0")]);
        assert!(given.state.is_active("var_synonyms"));
        assert_eq!(given.command_line, ["--var_synonyms"]);
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

    /// `--pick_order` is a `@LIST_FLAGS` entry read with Perl's `split`
    /// (`Config.pm:331`, `:761`), so a trailing comma names no criterion, with
    /// the value attached by `=` or given as the next token, and an empty value
    /// is the empty order, which `_is_flag_active` reads as inactive (`:902-906`)
    /// so the value check never meets it; a middle empty field is a value
    /// `%VALID` lacks (`:946`) and stays refused.
    #[test]
    fn pick_order_drops_a_trailing_empty_field_as_perl_does() {
        let canonical = Some(&["canonical".to_string()][..]);
        assert_eq!(
            config(&["vep", "--pick_order", "canonical,"])
                .pick_order
                .as_deref(),
            canonical
        );
        assert_eq!(
            config(&["vep", "--pick_order=canonical,"])
                .pick_order
                .as_deref(),
            canonical
        );
        assert_eq!(
            config(&["vep", "--pick_order", ""]).pick_order.as_deref(),
            Some(&[][..])
        );
        assert_eq!(
            config(&["vep", "--pick_order", ","]).pick_order.as_deref(),
            Some(&[][..])
        );
        let ParseFailure::Usage(error) =
            Config::from_argv(["vep", "--pick_order", "canonical,,rank"]).unwrap_err()
        else {
            panic!("a middle empty field is refused by the parser");
        };
        assert_eq!(error.kind(), clap::error::ErrorKind::InvalidValue);
        assert!(
            error.to_string().starts_with(
                "error: a value is required for '--pick_order <PICK_ORDER>' but none was supplied"
            ),
            "{error}"
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
            "ERROR: [cli] --af_exac was removed in Ensembl VEP release 107: use --af_gnomade"
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

    /// `transcript_filter` is read from the state after the option sets, which
    /// the `exclude_predicted` set appends to (Config.pm:566-571): the flag is
    /// in `@ALLOW_MULTIPLE` (Config.pm:340-346), so the set's filter is pushed
    /// after the command line's (Config.pm:870-871) wherever the flag stands.
    #[test]
    fn transcript_filter_follows_the_state_the_exclude_predicted_set_appends_to() {
        let predicted = "not stable_id match ^X._";
        assert_eq!(
            config(&["vep", "--exclude_predicted"]).transcript_filter,
            [predicted]
        );
        let biotype = "biotype is protein_coding";
        for argv in [
            ["vep", "--transcript_filter", biotype, "--exclude_predicted"],
            ["vep", "--exclude_predicted", "--transcript_filter", biotype],
        ] {
            assert_eq!(config(&argv).transcript_filter, [biotype, predicted]);
        }
        assert_eq!(
            config(&["vep", "--transcript_filter", biotype]).transcript_filter,
            [biotype]
        );
        assert!(config(&["vep"]).transcript_filter.is_empty());
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

    /// Ensembl VEP splits `--distance` on `,` with Perl's `split` (`Config.pm` 334,
    /// 761), which drops a trailing empty field, and reads the pair with `($u, $d) =
    /// @$distances; $d ||= $u` (`BaseRunner.pm` 505-506): one value serves both
    /// directions, `5000,` is that one value, `5000,2000` is asymmetric, the
    /// downstream part `0` is the upstream value (`2000,0` is `2000,2000`), an
    /// upstream 0 stands (`0,2000`), and `00` or `+0` is a true string in Perl, kept
    /// and numified to a downstream 0.
    #[test]
    fn distance_pair_with_zero_downstream_uses_upstream() {
        assert_eq!(parse_distance("5000"), Ok((5000, 5000)));
        assert_eq!(parse_distance("5000,"), Ok((5000, 5000)));
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

    /// A part that is not an integer, an empty part other than a trailing one, or a
    /// third part is refused with the contract line `ERROR: [config] <message>`
    /// naming the value as written. Ensembl VEP numifies a non-numeric part to 0
    /// under an `Argument "x" isn't numeric` warning, dropping the up- or downstream
    /// rows that part governs, ignores a third part, and reads ` 0 ` as a set
    /// downstream value (numified to 0) where a trimmed `0` would read as the
    /// upstream value, so a padded part is refused as well. A negative value in its
    /// own token (`--distance -1`) never reaches here: clap reads `-1` as an unknown
    /// flag and exits 2; `--distance=-1` is refused here.
    #[test]
    fn distance_non_numeric_is_an_error() {
        for value in [
            "x",
            "5000,x",
            "x,5000",
            "",
            ",",
            ",5000",
            "5000,,2000",
            "1.5",
            "1,2,3",
            " 2000",
            "2000, 0",
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

    /// `root_dir => dir_cache || dir` (AnnotationSourceAdaptor.pm:127): an empty
    /// or `0` `--dir_cache` is Perl-false and the root is `--dir`, or the
    /// directory `--cache` names (Config.pm:717's `||=` assigns over it too).
    #[test]
    fn test_dir_cache_empty_or_zero_falls_back_to_dir() {
        for value in ["", "0"] {
            assert_eq!(
                config(&["vep", "--dir", "/data/vep", "--dir_cache", value]).dir_cache,
                "/data/vep",
                "{value:?}"
            );
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().to_str().unwrap();
        assert_eq!(
            config(&["vep", "--cache", path, "--dir_cache", ""]).dir_cache,
            path
        );
    }

    /// The default `--dir` is `$ENV{HOME}.'/.vep'` (Config.pm:294): `/.vep`
    /// when `HOME` is unset, the home directory's `.vep` otherwise.
    #[test]
    fn test_dir_default_without_home_is_root_vep() {
        let expand = crate::argv::expand_tilde_text;
        assert_eq!(expand("~/.vep", ""), "/.vep");
        assert_eq!(expand("~/.vep", "/home/u"), "/home/u/.vep");
        assert_eq!(expand("~", ""), "/");
        assert_eq!(expand("~", "/home/u"), "/home/u");
        assert_eq!(expand("~x/.vep", ""), "~x/.vep");
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
        for component in [
            "cli",
            "input",
            "output",
            "cache",
            "fasta",
            "filter",
            "check_svs",
            "config",
            "run",
        ] {
            assert!(names.contains(component), "{component}");
        }
    }

    /// `--shift_hgvs 0` sets `shift_3prime` to 0 whatever the command line passes
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

    /// `--ga4gh_vrs` with `--vcf` is the pair Ensembl VEP refuses (Config.pm:629,
    /// `ERROR: Can't use --ga4gh_vrs and --vcf together`); with the other formats,
    /// and `--spdi` and `--hgvsg_use_accession` anywhere, the flags are plain
    /// switches that default off.
    #[test]
    fn ga4gh_vrs_with_vcf_is_an_error() {
        assert_eq!(
            refusal(&["vep", "--ga4gh_vrs", "--vcf"]),
            "ERROR: [cli] Can't use --ga4gh_vrs and --vcf together"
        );
        assert_eq!(
            refusal(&["vep", "--vcf", "--spdi", "--ga4gh_vrs"]),
            "ERROR: [cli] Can't use --ga4gh_vrs and --vcf together"
        );
        let c = config(&["vep", "--ga4gh_vrs", "--tab"]);
        assert!(c.ga4gh_vrs && !c.spdi);
        let c = config(&[
            "vep",
            "--ga4gh_vrs",
            "--json",
            "--spdi",
            "--hgvsg_use_accession",
        ]);
        assert!(c.ga4gh_vrs && c.spdi && c.hgvsg_use_accession);
        assert!(Config::from_argv(["vep", "--spdi", "--vcf"]).is_ok());
        let plain = config(&["vep"]);
        assert!(!plain.spdi && !plain.ga4gh_vrs && !plain.hgvsg_use_accession);
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
    /// (Runner.pm:907, whose `||` passes an empty string and `0` over as
    /// Perl-false), a leading `~` expanded as `glob` expands it (`:908`).
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
        for value in ["", "0"] {
            assert_eq!(
                config(&["vep", "--dir", "/data/vep", "--dir_plugins", value]).dir_plugins,
                "/data/vep/Plugins",
                "{value:?}"
            );
        }
        assert_eq!(
            config(&["vep", "--dir_plugins", "~/plugins"]).dir_plugins,
            format!("{home}/plugins")
        );
        assert!(!config(&["vep"]).safe);
        assert!(config(&["vep", "--safe"]).safe);
    }

    /// The warnings file is `--warning_file`, else `<output_file>_warnings.txt`
    /// with the output file as the user spelled it, `STDOUT` included
    /// (BaseVEP.pm:890 reads the `output_file` parameter as given), or
    /// `vep_warnings.txt` for an output named `0`, which Perl reads as false
    /// (`($self->param('output_file') || 'vep')`, the same line); a
    /// `--warning_file` that is false in Perl, the empty string or `0`, is no
    /// file and the default stands (`$file ||= ...`, the same line).
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
        assert_eq!(config(&["vep", "-o", "0"]).warning_file, "vep_warnings.txt");
        assert_eq!(
            config(&["vep", "-o", "x", "--warning_file", "w.txt"]).warning_file,
            "w.txt"
        );
        for false_value in ["", "0"] {
            assert_eq!(
                config(&["vep", "-o", "x", "--warning_file", false_value]).warning_file,
                "x_warnings.txt",
                "{false_value:?}"
            );
        }
        assert_eq!(
            config(&["vep", "--skipped_variants_file", "s.txt"])
                .skipped_variants_file
                .as_deref(),
            Some("s.txt")
        );
    }

    /// `--individual` sets `allow_non_variant` through its option set
    /// (Config.pm:413-417) and refuses `--minimal` (:626, the message of :964),
    /// in any order and under `--safe` not at all (:961); its value is a
    /// comma-separated list (:328) kept as given, `all` included.
    /// `--process_ref_homs` (:143) is a bare flag and implies nothing.
    #[test]
    fn individual_implies_allow_non_variant_and_refuses_minimal() {
        let c = config(&["vep", "--individual", "all"]);
        assert_eq!(
            c.individual.as_deref(),
            Some(["all".to_string()].as_slice())
        );
        assert!(c.allow_non_variant);
        assert!(!c.process_ref_homs);
        let c = config(&["vep", "--individual", "S1,S3", "--process_ref_homs"]);
        assert_eq!(
            c.individual.as_deref(),
            Some(["S1".to_string(), "S3".to_string()].as_slice())
        );
        assert!(c.allow_non_variant && c.process_ref_homs);
        let c = config(&["vep", "--process_ref_homs"]);
        assert!(c.process_ref_homs && !c.allow_non_variant && c.individual.is_none());
        assert_eq!(
            refusal(&["vep", "--individual", "all", "--minimal"]),
            "ERROR: [cli] Can't use --individual and --minimal together"
        );
        assert_eq!(
            refusal(&["vep", "--minimal", "--individual", "S1"]),
            "ERROR: [cli] Can't use --individual and --minimal together"
        );
        let safe = config(&["vep", "--safe", "--individual", "all", "--minimal"]);
        assert!(safe.minimal && safe.individual.is_some());
        assert_eq!(FlagState::implied_by("individual"), ["allow_non_variant"]);
    }

    /// `--individual_zyg` is a comma-separated list kept as given (Config.pm:329),
    /// sets `allow_non_variant` through its option set (:420-424) and is refused
    /// beside `--individual` (:626, the message of :964 naming the pair in the
    /// table's order), in any order and under `--safe` not at all (:961); it is
    /// compatible with `--minimal`, which the table lists against `--individual`
    /// alone.
    #[test]
    fn individual_zyg_refuses_individual() {
        let c = config(&["vep", "--individual_zyg", "all"]);
        assert_eq!(
            c.individual_zyg.as_deref(),
            Some(["all".to_string()].as_slice())
        );
        assert!(c.allow_non_variant && c.individual.is_none());
        let c = config(&["vep", "--individual_zyg", "S1,S3", "--minimal"]);
        assert_eq!(
            c.individual_zyg.as_deref(),
            Some(["S1".to_string(), "S3".to_string()].as_slice())
        );
        assert!(c.minimal);
        assert!(config(&["vep"]).individual_zyg.is_none());
        assert_eq!(
            refusal(&["vep", "--individual", "all", "--individual_zyg", "all"]),
            "ERROR: [cli] Can't use --individual and --individual_zyg together"
        );
        assert_eq!(
            refusal(&["vep", "--individual_zyg", "S1", "--individual", "S1"]),
            "ERROR: [cli] Can't use --individual and --individual_zyg together"
        );
        let safe = config(&[
            "vep",
            "--safe",
            "--individual",
            "all",
            "--individual_zyg",
            "all",
        ]);
        assert!(safe.individual.is_some() && safe.individual_zyg.is_some());
        assert_eq!(
            FlagState::implied_by("individual_zyg"),
            ["allow_non_variant"]
        );
    }
}
