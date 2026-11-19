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
use crate::input_shape::{CpuTopology, InputShape};
use crate::nearest::NearestType;
use crate::option_rules::{split_list_flag, FlagState, DEPRECATED, OWN_RULES};

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
    /// Whether the refusal applies to an active flag alone. Ensembl VEP refuses
    /// a deprecated flag only when its value is true (`Config.pm:977, 902-906`),
    /// so `--convert 0` proceeds; every other refused flag is refused when given.
    fn is_refused_when_active_alone(&self) -> bool {
        DEPRECATED.iter().any(|(flag, _)| *flag == self.flag)
    }

    /// The help text `--help` and the CLI reference print for the flag.
    pub fn help(&self) -> String {
        match self.treatment {
            Treatment::Refused { message, .. } if self.is_refused_when_active_alone() => {
                format!("Refused with an error unless its value is 0 or empty: {message}")
            }
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

/// `a` then `b` as bytes, in an array the caller sizes to their total length.
const fn joined<const N: usize>(a: &str, b: &str) -> [u8; N] {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    assert!(N == a.len() + b.len(), "the array holds both strings");
    let mut bytes = [0u8; N];
    let mut i = 0;
    while i < a.len() {
        bytes[i] = a[i];
        i += 1;
    }
    let mut j = 0;
    while j < b.len() {
        bytes[a.len() + j] = b[j];
        j += 1;
    }
    bytes
}

/// A refusal ending in [`GFF_GTF_RECIPE`], joined at compile time: the error
/// table's messages are `&'static str`.
macro_rules! with_recipe {
    ($prefix:expr) => {{
        const PREFIX: &str = $prefix;
        const BYTES: [u8; PREFIX.len() + GFF_GTF_RECIPE.len()] = joined(PREFIX, GFF_GTF_RECIPE);
        match std::str::from_utf8(&BYTES) {
            Ok(message) => message,
            Err(_) => panic!("two UTF-8 strings joined"),
        }
    }};
}

const GFF_MESSAGE: &str = with_recipe!(
    "--gff is not supported: Ensembl VEP reads a GFF3 as its transcript set through --custom file=<gff>,format=gff; "
);
const GTF_MESSAGE: &str = with_recipe!(
    "--gtf is not supported: Ensembl VEP reads a GTF as its transcript set through --custom file=<gtf>,format=gtf; "
);
const BAM_MESSAGE: &str = with_recipe!(
    "--bam is not supported: Ensembl VEP applies a BAM's alignments to RefSeq transcript models at run time; the RefSeq cache Ensembl distributes already carries those edits, which vep-rs reads from a converted homo_sapiens_refseq cache under --refseq; a GFF or GTF transcript set of your own is read without the BAM's edits; "
);

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
        "Config.pm:636-638 lists it as deprecated and :977-985 refuses it with this text when its value is true; a value Perl reads as false, `0` or the empty value, leaves the flag inactive (:902-906) and the run proceeds, there and here. Behind such a value, a true one from a lower source is refused there (:1101, 1158) and not read here",
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
        "gff",
        GFF_MESSAGE,
        "Config.pm:260 (gff=s); the option set at :531-536 sets --custom file=%gff%,format=gff; AnnotationSource/File/GFF.pm:71-116 lists the feature types read",
        "a cache built from the GFF3 by vep-cache-builder --gff3, read with --json_cache",
    ),
    refused(
        "gtf",
        GTF_MESSAGE,
        "Config.pm:261 (gtf=s); the option set at :538-543 sets --custom file=%gtf%,format=gtf; AnnotationSource/File/GTF.pm:78-85 lists the feature types read",
        "a cache built from the GTF by vep-cache-builder --gtf, read with --json_cache",
    ),
    refused(
        "bam",
        BAM_MESSAGE,
        "Config.pm:253 (bam=s); the option set at :573-579 switches on use_transcript_ref and bam_edited; AnnotationType/Transcript.pm:421-631 applies the alignments to each transcript as SeqEdits; CacheDir.pm:412-416 reads a BAM-edited cache's info",
        "--refseq on a converted homo_sapiens_refseq cache, which carries the edits",
    ),
    noop(
        "condel",
        "Ensembl VEP declares it and reads nothing from it",
        "Config.pm:206 (condel=s), read by no module",
    ),
    noop(
        "debug",
        "Ensembl VEP reads it only to print the plugin directories and modules it loads; the annotation output is unchanged",
        "Config.pm:275; Runner.pm:920, 941",
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
        "no module of Ensembl VEP reads it, and vep-rs writes no temporary file",
        "Config.pm:259 (tmpdir=s), default /tmp at :304, read by no module",
    ),
    noop(
        "offline",
        "vep-rs always runs offline and never contacts a database",
        "Config.pm:245; the option set at :464-469 switches on --cache and off --database",
    ),
];

/// The value `--cache` carries when given without a directory
/// ([`Args::cache`]'s `default_missing_value`).
pub const CACHE_FLAG_ALONE: &str = "1";

/// How a GFF or GTF transcript model is annotated against: the text a
/// `--custom` string in `format=gff` or `format=gtf` transcript mode is
/// refused with ([`vep_builtins::custom::GFF_GTF_RECIPE`]).
pub const GFF_GTF_RECIPE: &str = vep_builtins::custom::GFF_GTF_RECIPE;

/// The keys a `--regulatory_gff` string takes (`AnnotationSourceAdaptor.pm:77-81`
/// at release/116.2), in the sorted order its error texts list them.
pub const REGULATORY_GFF_KEYS: [&str; 3] = ["emars", "file", "motifs"];

/// The files `--regulatory_gff` names (`AnnotationSourceAdaptor.pm:227-262`):
/// `key=value` pairs joined by commas, or a bare file name, which is `file`.
/// Each file is a bgzip-compressed GFF with its tabix index. A value that is
/// false in Perl (`0` or empty) names no file.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RegulatoryGff {
    /// `file=`: the regulatory features.
    pub file: Option<String>,
    /// `emars=`: regulatory features with an identifier space of their own,
    /// read as `file=` is.
    pub emars: Option<String>,
    /// `motifs=`: the motif features.
    pub motifs: Option<String>,
}

impl RegulatoryGff {
    /// Parses the flag's value with Ensembl VEP's texts (`AnnotationSourceAdaptor.pm:237-262`):
    /// a string without `=` is the `file`; otherwise every comma-separated part,
    /// as Perl's `split` returns them ([`vep_core::sv_routing::perl_split`]: a
    /// trailing comma is no part), is `key=value` with a non-empty value, every
    /// key is one of [`REGULATORY_GFF_KEYS`], and at least one value is true in
    /// Perl.
    pub fn parse(spec: &str) -> Result<Self, CliError> {
        if !spec.contains('=') {
            return Ok(Self {
                file: Some(spec.to_string()),
                ..Self::default()
            });
        }
        let mut given: std::collections::BTreeMap<&str, &str> = Default::default();
        for param in vep_core::sv_routing::perl_split(spec, &[',']) {
            match param.split_once('=') {
                Some((key, value)) if !value.is_empty() => {
                    given.insert(key, value);
                }
                _ => {
                    return Err(CliError::cli(format!(
                        "Failed to parse --regulatory_gff parameter \"{param}\"; expected key=value"
                    )))
                }
            }
        }
        let invalid: Vec<&str> = given
            .keys()
            .copied()
            .filter(|key| !REGULATORY_GFF_KEYS.contains(key))
            .collect();
        if !invalid.is_empty() {
            return Err(CliError::cli(format!(
                "Unsupported --regulatory_gff key(s): {}\nSupported keys: {}",
                invalid.join(", "),
                REGULATORY_GFF_KEYS.join(", ")
            )));
        }
        let value = |key: &str| {
            given
                .get(key)
                .copied()
                .filter(|value| crate::argv::truthy(value))
                .map(str::to_string)
        };
        let parsed = Self {
            file: value("file"),
            emars: value("emars"),
            motifs: value("motifs"),
        };
        if parsed.sources().next().is_none() {
            return Err(CliError::cli(format!(
                "--regulatory_gff requires at least one of: {}",
                REGULATORY_GFF_KEYS.join(", ")
            )));
        }
        Ok(parsed)
    }

    /// The files named, each with its key, in the order Ensembl VEP builds its
    /// sources (`AnnotationSourceAdaptor.pm:284-310`): `file`, `emars`, `motifs`.
    pub fn sources(&self) -> impl Iterator<Item = (&'static str, &str)> {
        [
            ("file", self.file.as_deref()),
            ("emars", self.emars.as_deref()),
            ("motifs", self.motifs.as_deref()),
        ]
        .into_iter()
        .filter_map(|(key, file)| file.map(|file| (key, file)))
    }
}

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
        "parquet",
        "the Parquet writer: input out of chromosome and position order without --no_check_variants_order, or an output directory holding files vep-rs did not write",
    ),
    (
        "plugin",
        "a `--plugin` naming an Ensembl VEP plugin vep-rs does not provide, or a Perl module (a `.pm` file by path or under the plugins directories), with or without --safe; under --safe, a plugin that cannot be loaded or a plugins directory that does not exist",
    ),
    (
        "cache",
        "the cache: nothing resolves at the configured path, its layout or metadata is unusable, a transcript or regulatory shard cannot be read (a transcript an HGVS input names included), or a regulatory chromosome directory holds files but no shard, its assembly disagrees with --assembly or with its own directory name, or a variation file cannot be opened or read",
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
        "custom",
        "a `--custom` source: its string is malformed (Ensembl VEP's text), it names what vep-rs does not read (a remote URL, or a GFF or GTF source in transcript mode), or its file cannot be opened or read",
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
    let mut state = FlagState::from_matches(&command, &matches);
    // Getopt::Long stores an `=i` value as its number, and the option rules read
    // that for truth (Config.pm:902-906), so `--shift_3prime 00` is not active.
    for (flag, number) in [
        ("shift_3prime", args.shift_3prime),
        ("shift_genomic", args.shift_genomic),
    ] {
        if let Some(number) = number {
            state.set(flag, &number.to_string());
        }
    }
    Ok(ParsedArgs {
        args,
        state,
        command_line,
    })
}

/// Fully resolved configuration for a VEP run.
#[derive(Debug, Clone)]
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
    pub compress_output: Option<String>,
    /// `--parquet_partition`: the column the Parquet parts are laid out in
    /// directories by, `chrom` alone; `None` lays every part at the root.
    pub parquet_partition: Option<String>,
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
    /// The plugins directories: `--dir_plugins`, else `<dir>/Plugins`
    /// (Runner.pm:907), an empty or `0` value counting as absent, expanded as
    /// `glob` expands the value (`:908`) and split by whether each match is a
    /// directory (`:910-923`).
    pub dir_plugins: PluginDirs,
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
    /// The annotator pool's threads, as `--fork` was given; [`Fork::resolve`]
    /// turns `0` into a count once the input's shape and its inflater count are
    /// known.
    pub fork: Fork,
    /// BGZF decompression worker threads for a BGZF input in any format, as the
    /// flag was given (a plain gzip input decompresses on one thread);
    /// [`resolve_decompression_threads`] turns it into a count once the input's
    /// shape is known.
    pub decompression_threads: DecompressionThreads,
    /// BGZF deflate worker threads for `--compress_output bgzip`
    /// (`--compression_threads`, resolved; one is the writing thread itself).
    pub compression_threads: usize,
    /// Threads the cache shards parse on (`--parse_threads`); `None` is the
    /// pool's built-in width.
    pub parse_threads: Option<usize>,
    /// How far ahead of the reader the prewarm helper loads chromosomes
    /// (`--prewarm_ahead`).
    pub prewarm_ahead: PrewarmAhead,
    /// Whether the reader thread gets a physical core of its own on a host
    /// with hardware threads (`--reader_core`); [`reserve_reader_core`] decides
    /// once the input's record width is known.
    pub reader_core: ReaderCore,
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
    /// `--hgvsg_use_accession`: HGVSg names the sequence by the RefSeq synonym,
    /// an `NC_` name before any other, of the input name or of the sequence it
    /// resolves to, when the synonyms table has one.
    pub hgvsg_use_accession: bool,
    /// `--spdi`: the SPDI column.
    pub spdi: bool,
    /// The SPDI names its sequence as the cache names it on every row, as a
    /// database slice names itself (ensembl-variation VariationFeature.pm
    /// 2151): set by `vep recode`, whose reference reads a database. Off, a
    /// row that overlaps no transcript takes the FASTA's name, as Ensembl
    /// VEP's offline run does.
    pub spdi_on_the_cache_name: bool,
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
    /// `--cell_type` as given, when Ensembl VEP's comma split of it names
    /// anything (Config.pm:761; a value of commas alone is an empty list, which
    /// `_is_flag_active` reads as unset, :905); `None` otherwise. The names are
    /// read from it against the cache's list once the cache is resolved
    /// ([`split_cell_types`]), since an epigenome name may hold a comma.
    /// Switches `regulatory` on through the option set (`:427-432`).
    pub cell_type: Option<String>,
    /// `--regulatory_gff` (`Config.pm:210`): the regulatory and motif GFF files
    /// read in place of the cache's regulatory data; `None` without the flag or
    /// with a value that is false in Perl. The flag's option set switches
    /// `regulatory` on (`:434-439`).
    pub regulatory_gff: Option<RegulatoryGff>,
    /// `--extended_promoters` (`Config.pm:211`): a promoter read from
    /// `--regulatory_gff` is widened to its `extended_start` and `extended_end`.
    /// Only the GFF source reads it, so without `--regulatory_gff` it changes
    /// nothing.
    pub extended_promoters: bool,
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
    /// The HTML statistics report: the flag as given, or the default when
    /// neither `--stats_text` nor `--stats_html` names a report
    /// (`Runner.pm:787-790`); `--no_stats` overrides either.
    pub stats_html: bool,

    pub verbose: bool,
    pub quiet: bool,
    /// `--no_progress`: the progress line is not written; the status lines, the
    /// summary line and the output are as without it.
    pub no_progress: bool,
    /// The `INFO:` lines the resolution owes the run, in Ensembl VEP's order
    /// (`Runner.pm:117-120` before `post_setup_checks`), printed by the runner
    /// unless `--quiet`.
    pub notices: Vec<&'static str>,

    pub plugin: Vec<String>,
    pub custom: Vec<String>,
    /// `--custom_multi_allelic`: a `--custom` VCF source keeps a comma-separated
    /// INFO value whole on every allele (`File/VCF.pm:199`).
    pub custom_multi_allelic: bool,
    /// `--custom_suppress_filter`: a `--custom` VCF source adds its FILTER column
    /// only when its `fields` name it (`File/VCF.pm:219`).
    pub custom_suppress_filter: bool,

    pub refseq: bool,
    pub merged: bool,
    /// `--all_refseq` (`Config.pm:98`): every model of a RefSeq or merged
    /// cache is annotated, not the RefSeq accessions alone
    /// (`AnnotationType/Transcript.pm:237`); nothing changes on an Ensembl
    /// cache, where the RefSeq clause never runs.
    pub all_refseq: bool,
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
    /// `--no_check_variants_order` (`Config.pm:83`): the Parquet writer skips its
    /// order check and relaxes its layout; every text format ignores the flag.
    pub no_check_variants_order: bool,
    /// `--chr` as given, split on commas (`Config.pm:759-761`); `None` without the
    /// flag. The names it admits are read by `validate::ChromosomeList`.
    pub chr: Option<Vec<String>>,
    /// `--check_ref`: the reference allele is checked against the FASTA.
    pub check_ref: bool,
    /// `--lookup_ref`: the reference allele is read from the FASTA.
    pub lookup_ref: bool,
    /// `--use_transcript_ref` (`Config.pm:254`): each row's reference allele is
    /// the transcript's, and `GIVEN_REF` and `USED_REF` are written; off again
    /// under `--use_given_ref` (`:588-590`), and switched on by a BAM-edited
    /// cache unless that flag was given (`CacheDir.pm:412-416`).
    pub use_transcript_ref: bool,
    /// `--use_given_ref` (`Config.pm:255`) was given: the input's reference allele
    /// stands, and a BAM-edited cache leaves `use_transcript_ref` off.
    pub use_given_ref: bool,
    /// The cache's transcript models were edited against its alignment file,
    /// set by the cache resolution for a cache whose metadata names one, as
    /// `CacheDir.pm:415` sets it: the `BAM_EDIT` column is written. Never set by
    /// a flag: `--bam`, whose option set switches it on in Ensembl VEP
    /// (`Config.pm:573-579`), is refused.
    pub bam_edited: bool,
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

/// The remote root `--phyloP` and `--phastCons` read their bigWigs from by
/// default (`Config.pm:314`) is refused: the sources are read from a local
/// copy of the same layout, or named through `--custom`.
fn check_ucsc_data_root(state: &FlagState) -> Result<(), CliError> {
    let shortcut_given = state.is_active("phyloP") || state.is_active("phastCons");
    let remote_root = state
        .scalar("ucsc_data_root")
        .is_some_and(vep_builtins::custom::grammar::is_remote);
    if shortcut_given && remote_root {
        return Err(CliError::new(
            "custom",
            "--ucsc_data_root must name a local directory laid out as \
             <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw; the same source is \
             --custom file=<local>.bw,format=bigwig,type=exact",
        ));
    }
    Ok(())
}

/// The `--custom` strings the state holds after the option sets, the ones the
/// `phyloP` set appended with their short name spelt `phyloP<N>way`: Ensembl
/// VEP's set writes `phlyoP<N>way` (`Config.pm:555`), a misspelling whose
/// column vep-rs does not reproduce.
fn custom_sources(state: &FlagState) -> Vec<String> {
    // An unset flag leaves its `%name%` reference in the set's value as written
    // (`FlagState::expand`, `Config.pm:867`), so the string the set appended is
    // rebuilt here byte for byte and matched whole.
    let root = state.scalar("ucsc_data_root").unwrap_or("%ucsc_data_root%");
    let assembly = state.scalar("ucsc_assembly").unwrap_or("%ucsc_assembly%");
    let corrections: Vec<(String, String)> = state
        .list("phyloP")
        .iter()
        .map(|way| {
            let file = format!("{root}{assembly}/phyloP{way}way/{assembly}.phyloP{way}way.bw");
            (
                format!("{file},phlyoP{way}way,bigwig,exact"),
                format!("{file},phyloP{way}way,bigwig,exact"),
            )
        })
        .collect();
    state
        .list("custom")
        .into_iter()
        .map(|source| {
            corrections
                .iter()
                .find(|(as_set, _)| *as_set == source)
                .map_or(source, |(_, corrected)| corrected.clone())
        })
        .collect()
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
            let refusal_applies = if entry.is_refused_when_active_alone() {
                state.is_active(entry.flag)
            } else {
                state.is_given(entry.flag)
            };
            if refusal_applies {
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
        check_ucsc_data_root(&state)?;
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
        // Parsed apart from `resolve`, which cannot fail: a value false in Perl
        // names no file, and any other value must parse.
        let regulatory_gff = state
            .scalar("regulatory_gff")
            .filter(|value| crate::argv::truthy(value))
            .map(RegulatoryGff::parse)
            .transpose()?;
        let mut config = Config::resolve(args, &state, distance, delimiter);
        config.regulatory_gff = regulatory_gff;
        config.notices = notices;
        config.command_line = command_line;
        Ok(config)
    }

    /// The root the statistics report files are named from
    /// (`BaseRunner.pm:254`): `--stats_file` when its value is true in Perl, else
    /// `<output>_summary` on the output file as the command line spelled it, so
    /// `-o stdout` names `stdout_summary` in the working directory.
    pub fn stats_file_root(&self) -> String {
        crate::stats::stats_file_root(self.stats_file.as_deref(), &self.output_file_as_spelled)
    }

    /// Whether the input is standard input: `-i` names `STDIN` in any letter
    /// case (`uc($file) eq 'STDIN'`, Parser.pm:235). `input_file` keeps the
    /// spelling given, which the statistics report's `Input file` row
    /// (Stats.pm:652) and the skipped-variants file (BaseRunner.pm:371-372)
    /// print.
    pub fn reads_stdin(&self) -> bool {
        self.input_file.eq_ignore_ascii_case("STDIN")
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
        let plugins_pattern = args
            .dir_plugins
            .filter(|d| crate::argv::truthy(d))
            .unwrap_or_else(|| format!("{dir}/Plugins"));
        let dir_plugins = PluginDirs::expand(&plugins_pattern);

        // `cache_version || db_version || software_version` (CacheDir.pm:323):
        // a `0` is Perl-false and falls through to the release.
        let cache_version = args
            .cache_version
            .filter(|version| *version != 0)
            .unwrap_or(vep_core::VEP_VERSION);

        let fields = args.fields.as_deref().map(split_list_flag);

        // `shift_hgvs` is a setting only as 0 or 1 (BaseRunner.pm:498), and 0 sets
        // `shift_3prime` to 0 whatever the command line passes (Runner.pm:792-795).
        // Every reader of `shift_3prime` and `shift_genomic` tests the value for
        // truth (AnnotationType/Transcript.pm:162, OutputFactory.pm:1349, 1992), so
        // any integer but 0 is on.
        let shift_hgvs = args.shift_hgvs.as_deref() != Some("0");
        let is_on = |value: Option<i64>| value.is_some_and(|number| number != 0);
        let shift_3prime = shift_hgvs && is_on(args.shift_3prime);

        Config {
            input_file: args.input_file,
            input_data,
            format: state.scalar("format").unwrap_or("guess").to_string(),
            delimiter,
            output_file,
            output_file_as_spelled,
            output_format,
            compress_output: args.compress_output,
            parquet_partition: args.parquet_partition,
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
            shift_genomic: is_on(args.shift_genomic),
            shift_length: args.shift_length,
            species: args.species,
            assembly: state.scalar("assembly").map(str::to_string),
            buffer_size: args.buffer_size,
            max_sv_size: u64::try_from(args.max_sv_size).ok(),
            fork: Fork::from_raw(args.fork),
            // clap validated the text with `decompression_threads_text`.
            decompression_threads: match args.decompression_threads.as_deref() {
                Some(text) => text.parse().unwrap_or_default(),
                None => DECOMPRESSION_THREADS_DEFAULT,
            },
            compression_threads: resolve_compression_threads(args.compression_threads),
            parse_threads: resolve_parse_threads(args.parse_threads),
            // clap validated the text with `prewarm_ahead_text`.
            prewarm_ahead: args.prewarm_ahead.parse().unwrap_or_default(),
            // clap validated the text with `reader_core_text`.
            reader_core: args.reader_core.parse().unwrap_or_default(),
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
            spdi_on_the_cache_name: false,
            ga4gh_vrs: args.ga4gh_vrs,
            hgvsp_use_prediction: args.hgvsp_use_prediction,
            remove_hgvsp_version: args.remove_hgvsp_version,
            sift: state.scalar("sift").map(str::to_string),
            polyphen: state.scalar("polyphen").map(str::to_string),
            humdiv: state.scalar("polyphen_analysis") == Some("humdiv"),
            numbers: on("numbers") || vcf_output,
            domains: on("domains"),
            regulatory: on("regulatory"),
            cell_type: args.cell_type.filter(|_| on("cell_type")),
            regulatory_gff: None,
            extended_promoters: args.extended_promoters,
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
            stats_html: on("stats_html") || !args.stats_text,
            verbose: args.verbose,
            quiet: args.quiet,
            no_progress: args.no_progress,
            notices: Vec::new(),
            plugin: args.plugin,
            custom: custom_sources(state),
            custom_multi_allelic: args.custom_multi_allelic,
            custom_suppress_filter: args.custom_suppress_filter,
            refseq: args.refseq,
            merged: args.merged,
            all_refseq: args.all_refseq,
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
            no_check_variants_order: args.no_check_variants_order,
            chr: state.is_given("chr").then(|| state.list("chr")),
            check_ref: args.check_ref,
            lookup_ref: args.lookup_ref,
            use_transcript_ref: on("use_transcript_ref"),
            use_given_ref: args.use_given_ref,
            bam_edited: false,
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
pub const FORK_AUTO_CAP: usize = 32;

/// Resolve the --fork value. `0` means "auto": use the number of available
/// logical CPUs, capped at `FORK_AUTO_CAP`. Any non-zero value is returned
/// as-is so users can still pin the worker count explicitly.
///
/// The host's rule alone: a run that knows its input's shape resolves `0`
/// through [`Fork::resolve`], which lowers the count on an inflate-bound
/// input. Public so an embedding consumer of this crate can apply the same
/// auto-detect default without going through the full `Config::from_argv`
/// pipeline.
pub fn resolve_fork(raw: usize) -> usize {
    Fork::from_raw(raw).host_count()
}

/// The `--fork` value: `0` and an absent flag are `Auto`, resolved from the
/// host and the input's shape; a positive count is kept as given.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Fork {
    /// `0`: [`resolve_fork_auto`]'s rule.
    #[default]
    Auto,
    /// A positive count, uncapped.
    Fixed(usize),
}

impl Fork {
    /// The flag's value: `0` is `Auto`.
    pub fn from_raw(raw: usize) -> Self {
        if raw == 0 {
            Fork::Auto
        } else {
            Fork::Fixed(raw)
        }
    }

    /// The annotator count for a run over an input of `shape`, read by
    /// `inflaters` block-inflating threads (0 when the input is inflated on the
    /// reader's own thread or not at all), on this host.
    pub fn resolve(self, shape: &InputShape, inflaters: usize) -> usize {
        match self {
            Fork::Fixed(count) => count,
            Fork::Auto => resolve_fork_auto(shape, inflaters, fork_host_parallelism()),
        }
    }

    /// The annotator count without an input shape, the one an embedding
    /// consumer gets ([`resolve_fork`]).
    pub fn host_count(self) -> usize {
        self.resolve(&InputShape::Unknown, 0)
    }
}

/// The logical CPUs the `--fork` rule reads: 1 on a host that reports no count
/// (where the decompression rules assume [`UNREPORTED_HOST_CPUS`]).
fn fork_host_parallelism() -> usize {
    reported_parallelism().unwrap_or(1)
}

/// The `auto` rule for `--fork` on a host of `available_parallelism` logical
/// CPUs: the CPUs capped at [`FORK_AUTO_CAP`], except on an inflate-bound input
/// ([`is_inflate_bound`]), where the annotators take the CPUs the `inflaters`
/// leave, capped the same way and never fewer than half the capped CPU count
/// (nor than 1). The half floor is for the host whose inflate-bound pool is as
/// wide as its CPUs ([`INFLATE_BOUND_INFLATERS`] on a 32-CPU aarch64 host),
/// where the inflaters leave nothing and one annotator would carry the whole
/// file.
pub fn resolve_fork_auto(
    shape: &InputShape,
    inflaters: usize,
    available_parallelism: usize,
) -> usize {
    let host = available_parallelism.clamp(1, FORK_AUTO_CAP);
    if !is_inflate_bound(shape) {
        return host;
    }
    let beside_the_inflaters = available_parallelism.saturating_sub(inflaters);
    beside_the_inflaters.clamp((host / 2).max(1), FORK_AUTO_CAP)
}

/// An input whose records inflate to at least this many bytes is bound by
/// inflation rather than annotation; a wider inflater pool only pays there.
pub const INFLATE_BOUND_BYTES_PER_RECORD: u64 = 8_000;

/// Whether `shape` is an indexed input whose records inflate to at least
/// [`INFLATE_BOUND_BYTES_PER_RECORD`] bytes, the shape the `auto` rules of
/// `--decompression_threads` and `--fork` act on.
pub fn is_inflate_bound(shape: &InputShape) -> bool {
    matches!(
        shape,
        InputShape::Indexed {
            inflated_bytes_per_record,
            ..
        } if *inflated_bytes_per_record >= INFLATE_BOUND_BYTES_PER_RECORD
    )
}

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompressionThreads {
    /// `0`: the host's available parallelism, capped at
    /// [`DECOMPRESSION_THREADS_AUTO_CAP`].
    HostCapped,
    /// `auto`: the count the input's shape asks for
    /// ([`resolve_decompression_threads_auto`]).
    Auto,
    /// A positive count, lowered to the host's available parallelism.
    Fixed(usize),
}

/// The value an absent `--decompression_threads` takes: `auto` on an aarch64
/// host, where an inflate-bound input gains from the wider pool, and `0` on
/// any other, where [`INFLATE_BOUND_INFLATERS`] is the `0` rule's cap, so `auto`
/// would read the input's shape to resolve the same count.
#[cfg(target_arch = "aarch64")]
pub const DECOMPRESSION_THREADS_DEFAULT: DecompressionThreads = DecompressionThreads::Auto;

/// The value an absent `--decompression_threads` takes on a host that is not
/// aarch64: the `0` rule.
#[cfg(not(target_arch = "aarch64"))]
pub const DECOMPRESSION_THREADS_DEFAULT: DecompressionThreads = DecompressionThreads::HostCapped;

impl Default for DecompressionThreads {
    /// [`DECOMPRESSION_THREADS_DEFAULT`].
    fn default() -> Self {
        DECOMPRESSION_THREADS_DEFAULT
    }
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
    if is_inflate_bound(shape) {
        INFLATE_BOUND_INFLATERS.min(available_parallelism)
    } else {
        available_parallelism.min(DECOMPRESSION_THREADS_AUTO_CAP)
    }
}

/// Cap for the `0` rule of `--compression_threads`, the default: a worker
/// deflates a few hundred megabytes of output text per second, and a run's wall
/// time stops falling once about sixteen workers keep pace with the annotators'
/// output on a wide host, so the rule stops there rather than starting a deflate
/// thread per CPU.
pub const COMPRESSION_THREADS_AUTO_CAP: usize = 16;

/// Resolve the `--compression_threads` value: `0` is the host's available
/// logical CPUs capped at [`COMPRESSION_THREADS_AUTO_CAP`]; a positive count is
/// lowered to the host's CPUs when the host reports them (a worker beyond them
/// deflates nothing more), else kept as given. One is the single-threaded
/// encoder on the writing thread.
pub fn resolve_compression_threads(raw: usize) -> usize {
    resolve_compression_threads_on(raw, reported_parallelism())
}

/// [`resolve_compression_threads`] on a host whose logical CPU count is
/// `reported`, `None` when the host does not report one.
fn resolve_compression_threads_on(raw: usize, reported: Option<usize>) -> usize {
    match raw {
        0 => reported
            .unwrap_or(UNREPORTED_HOST_CPUS)
            .min(COMPRESSION_THREADS_AUTO_CAP),
        n => reported.map_or(n, |cpus| n.min(cpus)),
    }
}

/// The `--reader_core` value: whether the reader thread, the one that scans
/// the inflated input text, is given a physical core no other thread of the
/// run may use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReaderCore {
    /// `auto`: a core is reserved by [`reserve_reader_core`]'s rule once the
    /// input's records read as inflate-bound.
    #[default]
    Auto,
    /// `none`: no core is reserved.
    None,
    /// A logical CPU number: the reader runs there and that CPU's core is
    /// reserved for it, whatever the input.
    Fixed(usize),
}

impl FromStr for ReaderCore {
    type Err = String;

    /// Accepts exactly `auto`, `none` and a CPU number written in digits alone.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "auto" {
            return Ok(ReaderCore::Auto);
        }
        if s == "none" {
            return Ok(ReaderCore::None);
        }
        let digits_only = !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        match s.parse::<usize>() {
            Ok(cpu) if digits_only => Ok(ReaderCore::Fixed(cpu)),
            _ => Err(format!("'{s}' is not auto, none or a logical CPU number")),
        }
    }
}

impl fmt::Display for ReaderCore {
    /// The value as the flag spells it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReaderCore::Auto => f.write_str("auto"),
            ReaderCore::None => f.write_str("none"),
            ReaderCore::Fixed(cpu) => write!(f, "{cpu}"),
        }
    }
}

/// The `--reader_core` text, validated as [`ReaderCore`] parses it; the
/// command-line parser's value parser.
pub fn reader_core_text(text: &str) -> Result<String, String> {
    text.parse::<ReaderCore>().map(|_| text.to_string())
}

/// Physical cores a host needs before the `auto` rule reserves one for the
/// reader: on a smaller host the annotators and the inflaters would lose too
/// large a share of it.
pub const READER_CORE_MIN_PHYSICAL_CORES: usize = 8;

/// The threads outside the annotator pool and the inflater pool that stay busy
/// through a run (the reader, the block reader's frame thread and the writer),
/// counted when the run is tested for oversubscribing the physical cores.
const READER_CORE_OTHER_THREADS: usize = 3;

/// Records the reader scans before it decides its core: enough lines for the
/// inflated bytes per record to stand for the file (a header's bytes spread
/// over them move the figure by a few bytes), and few enough that a run of
/// seconds decides within its first hundredth.
pub const READER_CORE_DECISION_RECORDS: u64 = 4_096;

/// The CPUs a reserved reader core leaves to each side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReaderCoreReservation {
    /// The CPU the reader runs on.
    pub reader_cpu: usize,
    /// The CPUs every other thread of the run keeps: the allowed set less the
    /// reader's core.
    pub others: Vec<usize>,
}

/// The reader's core under `value`, on a host of `topology` whose process may
/// run on the `allowed` CPUs, for a run of `annotators` pool threads and
/// `inflaters` block-inflating threads over an input whose records inflate to
/// `inflated_bytes_per_record` bytes.
///
/// `auto` reserves the highest-numbered core whose every CPU is allowed when
/// the host has hardware threads (a core with two or more allowed CPUs), the
/// allowed CPUs make up at least [`READER_CORE_MIN_PHYSICAL_CORES`] whole
/// cores, the run's threads outnumber those cores (`annotators + inflaters`
/// plus the reader, frame and writer threads), and the records inflate to at
/// least [`INFLATE_BOUND_BYTES_PER_RECORD`] bytes, the width at which the
/// reader's scan of the inflated text sets the wall and every thread sharing
/// its core slows it. A fixed CPU puts the reader on that CPU and reserves its
/// core when the CPU and its core's other CPUs are allowed and leave another
/// CPU for the rest of the run. `none`, a host without a topology, and any
/// failed test reserve nothing.
pub fn reserve_reader_core(
    value: ReaderCore,
    topology: Option<&CpuTopology>,
    allowed: &[usize],
    annotators: usize,
    inflaters: usize,
    inflated_bytes_per_record: u64,
) -> Option<ReaderCoreReservation> {
    let topology = topology?;
    let whole_core = |core: &[usize]| core.iter().all(|cpu| allowed.contains(cpu));
    let reserve = |core: &[usize], reader_cpu: usize| {
        let others: Vec<usize> = allowed
            .iter()
            .copied()
            .filter(|cpu| !core.contains(cpu))
            .collect();
        (!others.is_empty()).then_some(ReaderCoreReservation { reader_cpu, others })
    };
    match value {
        ReaderCore::None => None,
        ReaderCore::Fixed(cpu) => {
            let core = topology.core_of(cpu)?;
            if !whole_core(core) {
                return None;
            }
            reserve(core, cpu)
        }
        ReaderCore::Auto => {
            if inflated_bytes_per_record < INFLATE_BOUND_BYTES_PER_RECORD {
                return None;
            }
            let whole_cores: Vec<&[usize]> = topology
                .cores()
                .iter()
                .map(Vec::as_slice)
                .filter(|core| whole_core(core))
                .collect();
            let hardware_threads = whole_cores.iter().any(|core| core.len() >= 2);
            if !hardware_threads
                || whole_cores.len() < READER_CORE_MIN_PHYSICAL_CORES
                || annotators + inflaters + READER_CORE_OTHER_THREADS <= whole_cores.len()
            {
                return None;
            }
            whole_cores
                .iter()
                .rev()
                .find(|core| core.len() >= 2)
                .and_then(|core| reserve(core, core[0]))
        }
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

/// Resolve the `--parse_threads` value: `0` is `None`, the pool's built-in
/// width (the host's available parallelism), else the thread count the cache
/// shards parse on. The flag's parser bounds the value, so it fits a `usize`.
pub fn resolve_parse_threads(raw: u64) -> Option<usize> {
    if raw == 0 {
        None
    } else {
        Some(raw as usize)
    }
}

/// The `--prewarm_ahead` value: how far ahead of the reader the prewarm helper
/// loads the input's chromosomes, in the input index's contig order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrewarmAhead {
    /// A count of chromosomes after the one the reader has reached; `0` warms
    /// the reached chromosome alone.
    Chromosomes(usize),
    /// `all`: every indexed contig the cache holds, queued when the reader
    /// starts, so the helper loads the whole file's chromosomes in index order
    /// however far the reader has got.
    All,
}

/// The chromosomes an absent `--prewarm_ahead` warms ahead of the reader: one.
pub const PREWARM_AHEAD_DEFAULT: PrewarmAhead = PrewarmAhead::Chromosomes(1);

impl Default for PrewarmAhead {
    /// [`PREWARM_AHEAD_DEFAULT`].
    fn default() -> Self {
        PREWARM_AHEAD_DEFAULT
    }
}

impl FromStr for PrewarmAhead {
    type Err = String;

    /// Accepts exactly `all` and a count written in digits alone.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "all" {
            return Ok(PrewarmAhead::All);
        }
        let digits_only = !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        match s.parse::<usize>() {
            Ok(count) if digits_only => Ok(PrewarmAhead::Chromosomes(count)),
            _ => Err(format!("'{s}' is not all or a whole number of chromosomes")),
        }
    }
}

impl fmt::Display for PrewarmAhead {
    /// The value as the flag spells it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PrewarmAhead::Chromosomes(count) => write!(f, "{count}"),
            PrewarmAhead::All => f.write_str("all"),
        }
    }
}

/// The `--prewarm_ahead` text, validated as [`PrewarmAhead`] parses it; the
/// command-line parser's value parser.
pub fn prewarm_ahead_text(text: &str) -> Result<String, String> {
    text.parse::<PrewarmAhead>().map(|_| text.to_string())
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

/// The epigenomes a `--cell_type` value names, each as its position in the
/// cache's `cell_types` list, in the order given, a name given twice listed
/// twice. The value is read by greedy longest match: at each position the
/// longest list entry that the rest of the value starts with, followed by a
/// comma or the end, is the next name, and the comma is consumed. Ensembl VEP
/// splits the value on every comma (Config.pm:761) and tests each piece against
/// the list it splits the cache's comma-joined `cell_types` line on
/// (CacheDir.pm:174), so a name holding a comma (77 of the 112 GRCh38 names,
/// `A549_(m,_58_y)` among them) can never be asked for there. A remainder of
/// commas alone is the trailing empty fields Perl's `split` drops, and names
/// nothing more; any other remainder that starts with no name is Ensembl VEP's
/// refusal (AnnotationType/RegFeat.pm:85-92), naming the remainder up to its
/// next comma, the piece Ensembl VEP's split would name, over the list.
pub fn split_cell_types(value: &str, available: &[String]) -> Result<Vec<usize>, CliError> {
    let mut names = Vec::new();
    let mut rest = value;
    while !rest.is_empty() {
        let trailing_empty_fields = rest.bytes().all(|b| b == b',');
        if trailing_empty_fields {
            break;
        }
        let longest = available
            .iter()
            .enumerate()
            .filter(|(_, name)| {
                rest.starts_with(name.as_str())
                    && matches!(rest.as_bytes().get(name.len()), None | Some(b','))
            })
            .max_by_key(|(_, name)| name.len());
        let Some((index, name)) = longest else {
            let unknown = rest.split_once(',').map_or(rest, |(piece, _)| piece);
            return Err(CliError::cli(format!(
                "Cell type \"{unknown}\" unavailable; available cell types are:\n{}",
                available.join(",")
            )));
        };
        names.push(index);
        rest = rest[name.len()..].strip_prefix(',').unwrap_or("");
    }
    Ok(names)
}

/// Expand a leading `~` to the user's home directory
/// ([`crate::argv::expand_tilde_text`] with `HOME`, empty when it is unset).
fn expand_tilde(path: &str) -> String {
    crate::argv::expand_tilde_text(path, &std::env::var("HOME").unwrap_or_default())
}

/// The plugins directories of a run, from the one value `get_all_Plugins` reads
/// (Runner.pm:907-923): every path `glob` expands it to ([`crate::glob::expand`],
/// `:908`), each a search directory when it is one and a warning otherwise.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PluginDirs {
    /// The value as given: `--dir_plugins`, else `<dir>/Plugins`.
    pub pattern: String,
    /// The matches that are directories, the last match first: `unshift @INC,
    /// $dir` over the matches in order (`:921`) searches the last one first.
    pub dirs: Vec<String>,
    /// The matches that are not directories, in match order, each `Plugins
    /// directory '<match>' not found.` when a plugin is named (`:911-916`).
    pub not_found: Vec<String>,
}

impl PluginDirs {
    /// `pattern` expanded under `HOME` (empty when it is unset) and split.
    pub fn expand(pattern: &str) -> Self {
        let home = std::env::var("HOME").unwrap_or_default();
        let mut dirs = Vec::new();
        let mut not_found = Vec::new();
        for path in crate::glob::expand(pattern, &home) {
            if std::path::Path::new(&path).is_dir() {
                dirs.push(path);
            } else {
                not_found.push(path);
            }
        }
        dirs.reverse();
        Self {
            pattern: pattern.to_string(),
            dirs,
            not_found,
        }
    }

    /// Where a plugin that resolved to no file was looked for, for the `Failed to
    /// load plugin` line: the directories searched; with none, the one match
    /// that is no directory, or the value as given when there is no such match
    /// or several.
    pub fn searched(&self) -> String {
        match (self.dirs.as_slice(), self.not_found.as_slice()) {
            ([], [missing]) => missing.clone(),
            ([], _) => self.pattern.clone(),
            (dirs, _) => dirs.join(", "),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::option_rules::OPTION_SETS;
    use clap::CommandFactory;

    fn config(argv: &[&str]) -> Config {
        Config::from_argv(argv).unwrap()
    }

    /// An indexed input shape whose records inflate to the given bytes.
    fn indexed(inflated_bytes_per_record: u64) -> InputShape {
        InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate: 1_000,
            file_bytes: 1 << 20,
        }
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
    /// explicit value passes through uncapped, through `resolve_fork` and
    /// through the configuration's `Fork` alike.
    #[test]
    fn auto_fork_is_capped_and_explicit_fork_is_not() {
        let auto = resolve_fork(0);
        assert!(
            (1..=FORK_AUTO_CAP).contains(&auto),
            "auto fork {auto} outside 1..={FORK_AUTO_CAP}"
        );
        assert_eq!(resolve_fork(1), 1);
        assert_eq!(resolve_fork(64), 64);
        assert_eq!(config(&["vep"]).fork, Fork::Auto);
        assert_eq!(config(&["vep", "--fork", "0"]).fork, Fork::Auto);
        assert_eq!(config(&["vep", "--fork", "64"]).fork, Fork::Fixed(64));
        assert_eq!(Fork::Auto.host_count(), auto);
        assert_eq!(Fork::Fixed(64).host_count(), 64);
        assert_eq!(
            Fork::Fixed(3).resolve(&indexed(INFLATE_BOUND_BYTES_PER_RECORD), 10),
            3
        );
        assert_eq!(Fork::Auto.resolve(&InputShape::Unknown, 10), auto);
    }

    /// An absent `--decompression_threads` is `auto` on an aarch64 host and the
    /// `0` rule on any other, where `auto` would resolve the same count on every
    /// shape: `INFLATE_BOUND_INFLATERS` is the `0` rule's cap there.
    #[test]
    fn decompression_threads_default_is_auto_on_aarch64_and_host_capped_elsewhere() {
        let default = config(&["vep"]).decompression_threads;
        assert_eq!(default, DecompressionThreads::default());
        assert_eq!(default, DECOMPRESSION_THREADS_DEFAULT);
        if cfg!(target_arch = "aarch64") {
            assert_eq!(default, DecompressionThreads::Auto);
            assert_eq!(INFLATE_BOUND_INFLATERS, 32);
        } else {
            assert_eq!(default, DecompressionThreads::HostCapped);
            assert_eq!(INFLATE_BOUND_INFLATERS, DECOMPRESSION_THREADS_AUTO_CAP);
            for shape in [
                InputShape::Unknown,
                indexed(INFLATE_BOUND_BYTES_PER_RECORD - 1),
                indexed(INFLATE_BOUND_BYTES_PER_RECORD),
            ] {
                for cpus in [1, 4, 10, 32, 64] {
                    assert_eq!(
                        resolve_decompression_threads_on(
                            DecompressionThreads::Auto,
                            &shape,
                            Some(cpus)
                        ),
                        resolve_decompression_threads_on(
                            DecompressionThreads::HostCapped,
                            &shape,
                            Some(cpus)
                        ),
                        "{shape:?} on {cpus} CPUs"
                    );
                }
            }
        }
        // The explicit forms keep their meaning on every architecture.
        assert_eq!(
            config(&["vep", "--decompression_threads", "0"]).decompression_threads,
            DecompressionThreads::HostCapped
        );
        assert_eq!(
            config(&["vep", "--decompression_threads", "auto"]).decompression_threads,
            DecompressionThreads::Auto
        );
    }

    /// `--fork 0` on an inflate-bound input gives the annotators the CPUs the
    /// inflaters leave, capped at `FORK_AUTO_CAP` and floored at half the capped
    /// CPU count; on any other shape, and whatever the inflaters, it is the
    /// host's rule.
    #[test]
    fn fork_auto_leaves_the_inflaters_their_cores() {
        let wide = indexed(INFLATE_BOUND_BYTES_PER_RECORD);
        let narrow = indexed(INFLATE_BOUND_BYTES_PER_RECORD - 1);
        // (shape, inflaters, cpus, annotators)
        for (shape, inflaters, cpus, annotators) in [
            // The inflaters' CPUs are left to them.
            (&wide, 10, 32, 22),
            (&wide, 10, 16, 8),
            (&wide, 4, 8, 4),
            (&wide, 2, 8, 6),
            // The cap holds where the inflaters leave more than it.
            (&wide, 10, 64, 32),
            (&wide, 32, 64, 32),
            (&wide, 0, 64, 32),
            // The floor: an inflater pool as wide as the host leaves nothing.
            (&wide, 32, 32, 16),
            (&wide, 4, 4, 2),
            (&wide, 2, 2, 1),
            (&wide, 1, 1, 1),
            (&wide, 40, 32, 16),
            // A narrow or unknown shape keeps the host's rule.
            (&narrow, 10, 32, 32),
            (&narrow, 32, 32, 32),
            (&InputShape::Unknown, 10, 32, 32),
            (&InputShape::Unknown, 0, 64, 32),
            (&InputShape::Unknown, 0, 1, 1),
        ] {
            assert_eq!(
                resolve_fork_auto(shape, inflaters, cpus),
                annotators,
                "{shape:?} with {inflaters} inflaters on {cpus} CPUs"
            );
        }
        // Never below 1, whatever the host reports.
        assert_eq!(resolve_fork_auto(&wide, 5, 0), 1);
        assert_eq!(resolve_fork_auto(&InputShape::Unknown, 0, 0), 1);
        // The inflate-bound test is the one `--decompression_threads auto` uses.
        assert!(is_inflate_bound(&wide));
        assert!(!is_inflate_bound(&narrow));
        assert!(!is_inflate_bound(&InputShape::Unknown));
    }

    /// `--parse_threads` and `--prewarm_ahead` reach the configuration: `0`
    /// parse threads is the pool's built-in width, a count pins it and the
    /// flag's bound refuses more; `--prewarm_ahead` takes a count or `all`, one
    /// chromosome without the flag.
    #[test]
    fn parse_pool_is_bounded() {
        assert_eq!(resolve_parse_threads(0), None);
        assert_eq!(resolve_parse_threads(1), Some(1));
        assert_eq!(resolve_parse_threads(1024), Some(1024));
        assert_eq!(config(&["vep"]).parse_threads, None);
        assert_eq!(config(&["vep", "--parse_threads", "0"]).parse_threads, None);
        assert_eq!(
            config(&["vep", "--parse_threads", "16"]).parse_threads,
            Some(16)
        );
        for text in ["1025", "-1", "auto", "1.5"] {
            match Config::from_argv(["vep", "--parse_threads", text]).unwrap_err() {
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

        assert_eq!(config(&["vep"]).prewarm_ahead, PrewarmAhead::Chromosomes(1));
        assert_eq!(config(&["vep"]).prewarm_ahead, PrewarmAhead::default());
        for (text, value) in [
            ("0", PrewarmAhead::Chromosomes(0)),
            ("1", PrewarmAhead::Chromosomes(1)),
            ("2", PrewarmAhead::Chromosomes(2)),
            ("all", PrewarmAhead::All),
        ] {
            assert_eq!(text.parse::<PrewarmAhead>().unwrap(), value, "{text}");
            assert_eq!(prewarm_ahead_text(text).unwrap(), text);
            assert_eq!(value.to_string(), text);
            assert_eq!(
                config(&["vep", "--prewarm_ahead", text]).prewarm_ahead,
                value,
                "{text}"
            );
        }
        for text in ["-1", "ALL", "", "+1", "1.5", "every"] {
            assert!(text.parse::<PrewarmAhead>().is_err(), "{text}");
            assert!(prewarm_ahead_text(text).is_err(), "{text}");
            match Config::from_argv(["vep", "--prewarm_ahead", text]).unwrap_err() {
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
            DecompressionThreads::default()
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

    /// `--reader_core` takes exactly `auto`, `none` and a CPU number in digits;
    /// everything else is clap's usage error, and the default is `auto`.
    #[test]
    fn reader_core_parses_its_three_forms() {
        for (text, value) in [
            ("auto", ReaderCore::Auto),
            ("none", ReaderCore::None),
            ("0", ReaderCore::Fixed(0)),
            ("31", ReaderCore::Fixed(31)),
        ] {
            assert_eq!(text.parse::<ReaderCore>().unwrap(), value, "{text}");
            assert_eq!(reader_core_text(text).unwrap(), text);
            assert_eq!(
                config(&["vep", "--reader_core", text]).reader_core,
                value,
                "{text}"
            );
            assert_eq!(value.to_string(), text);
        }
        assert_eq!(config(&["vep"]).reader_core, ReaderCore::Auto);
        for text in ["-1", "1.5", "AUTO", "None", "", "+1", "0x10", "auto1"] {
            assert!(text.parse::<ReaderCore>().is_err(), "{text}");
            assert!(reader_core_text(text).is_err(), "{text}");
            match Config::from_argv(["vep", "--reader_core", text]).unwrap_err() {
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

    /// The reader-core rule over its inputs, on a sixteen-core host of paired
    /// hardware threads (CPU `i` with `i + 16`) with every CPU allowed:
    /// `auto` reserves the last core, CPUs 15 and 31, for an inflate-bound input
    /// of a run that oversubscribes the cores (16 annotators, 10 inflaters), the
    /// reader on CPU 15 and the other 30 CPUs to the rest; the same host reserves
    /// nothing for an input one byte under the threshold, for a run of 1
    /// annotator and 10 inflaters (14 threads on 16 cores), under `none`, or
    /// without a topology. A host whose cores are single CPUs (no hardware
    /// threads), one of 7 paired cores, and one whose allowed set lacks a CPU of
    /// every core but 7 reserve nothing under `auto`; one whose allowed set
    /// lacks CPU 31 alone reserves the last whole core, CPUs 14 and 30. A fixed
    /// CPU puts the reader on that CPU and reserves its core whatever the input
    /// and the counts, and nothing when the CPU is outside the topology, outside
    /// the allowed set, or its core is the only one allowed.
    #[test]
    fn reader_core_reserved_on_oversubscribed_paired_cores_only() {
        let paired = |cores: usize| {
            CpuTopology::from_cores((0..cores).map(|core| vec![core, core + cores]).collect())
        };
        let single = CpuTopology::from_cores((0..32).map(|cpu| vec![cpu]).collect());
        let host = paired(16);
        let all: Vec<usize> = (0..32).collect();
        let wide = INFLATE_BOUND_BYTES_PER_RECORD;
        let reserve = |value, topology: Option<&CpuTopology>, allowed: &[usize], a, i, b| {
            reserve_reader_core(value, topology, allowed, a, i, b)
        };

        let reserved = reserve(ReaderCore::Auto, Some(&host), &all, 16, 10, wide).unwrap();
        assert_eq!(reserved.reader_cpu, 15);
        assert_eq!(
            reserved.others,
            all.iter()
                .copied()
                .filter(|&c| c != 15 && c != 31)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            reserve(ReaderCore::Auto, Some(&host), &all, 16, 10, wide - 1),
            None
        );
        assert_eq!(
            reserve(ReaderCore::Auto, Some(&host), &all, 1, 10, wide),
            None
        );
        assert!(reserve(ReaderCore::Auto, Some(&host), &all, 4, 10, wide).is_some());
        assert_eq!(
            reserve(ReaderCore::None, Some(&host), &all, 16, 10, wide),
            None
        );
        assert_eq!(reserve(ReaderCore::Auto, None, &all, 16, 10, wide), None);
        assert_eq!(
            reserve(ReaderCore::Auto, Some(&single), &all, 16, 10, wide),
            None
        );
        let small = paired(7);
        let small_all: Vec<usize> = (0..14).collect();
        assert_eq!(
            reserve(ReaderCore::Auto, Some(&small), &small_all, 16, 10, wide),
            None
        );
        // Seven whole cores: CPUs 0 to 6 with 16 to 22, the other cores missing a CPU.
        let seven_whole: Vec<usize> = (0..7).chain(16..23).chain(7..16).collect();
        assert_eq!(
            reserve(ReaderCore::Auto, Some(&host), &seven_whole, 16, 10, wide),
            None
        );
        let without_31: Vec<usize> = (0..31).collect();
        let reserved = reserve(ReaderCore::Auto, Some(&host), &without_31, 16, 10, wide).unwrap();
        assert_eq!(reserved.reader_cpu, 14);
        assert_eq!(reserved.others.len(), 29);
        assert!(!reserved.others.contains(&30));

        let fixed = reserve(ReaderCore::Fixed(3), Some(&host), &all, 1, 0, 0).unwrap();
        assert_eq!(fixed.reader_cpu, 3);
        assert!(!fixed.others.contains(&3) && !fixed.others.contains(&19));
        assert_eq!(fixed.others.len(), 30);
        let fixed = reserve(ReaderCore::Fixed(19), Some(&host), &all, 1, 0, 0).unwrap();
        assert_eq!(fixed.reader_cpu, 19);
        assert!(!fixed.others.contains(&3) && !fixed.others.contains(&19));
        assert_eq!(
            reserve(ReaderCore::Fixed(99), Some(&host), &all, 1, 0, 0),
            None
        );
        assert_eq!(
            reserve(ReaderCore::Fixed(31), Some(&host), &without_31, 1, 0, 0),
            None
        );
        assert_eq!(
            reserve(ReaderCore::Fixed(0), Some(&host), &[0, 16], 1, 0, 0),
            None
        );
        assert_eq!(reserve(ReaderCore::Fixed(0), None, &all, 1, 0, 0), None);
    }

    /// The `auto` rule over the three shapes and both sides of the threshold:
    /// an unknown shape and an input below the threshold take the `0` rule, an
    /// input at or above it takes the inflate-bound pool, and the host's
    /// parallelism caps both.
    #[test]
    fn auto_decompression_threads_follow_the_input_shape() {
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

    /// The `0` compression worker count is the host's CPUs within its cap; a
    /// pinned count is kept up to the host's logical CPUs and lowered to them
    /// above that, and kept whole on a host that reports no count.
    #[test]
    fn compression_threads_pinned_up_to_the_host_cpus() {
        let cpus = std::thread::available_parallelism().unwrap().get();
        let host_capped = resolve_compression_threads(0);
        assert_eq!(host_capped, cpus.min(COMPRESSION_THREADS_AUTO_CAP));
        assert!((1..=COMPRESSION_THREADS_AUTO_CAP).contains(&host_capped));
        assert_eq!(resolve_compression_threads(1), 1);
        assert_eq!(resolve_compression_threads(cpus), cpus);
        assert_eq!(resolve_compression_threads(cpus + 1), cpus);
        assert_eq!(resolve_compression_threads(usize::MAX), cpus);

        assert_eq!(resolve_compression_threads_on(1, None), 1);
        assert_eq!(resolve_compression_threads_on(64, None), 64);
        assert_eq!(
            resolve_compression_threads_on(0, None),
            UNREPORTED_HOST_CPUS
        );
        assert_eq!(resolve_compression_threads_on(64, Some(16)), 16);
        assert_eq!(resolve_compression_threads_on(0, Some(2)), 2);
        assert_eq!(
            resolve_compression_threads_on(0, Some(COMPRESSION_THREADS_AUTO_CAP * 4)),
            COMPRESSION_THREADS_AUTO_CAP
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
    /// `--regulatory` part writes, the regulatory and motif rows with every
    /// motif field, and names no exclusion: `--everything` switches on exactly
    /// Ensembl VEP's set (Config.pm:351-380), `cell_type` outside it on both
    /// engines, so `CELL_TYPE` comes from `--cell_type` beside the set.
    #[test]
    fn everything_help_names_no_exclusion() {
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
            "every motif field included",
            "CELL_TYPE is added by --cell_type",
        ] {
            assert!(help.contains(phrase), "{phrase:?} missing from {help}");
        }
        for word in [
            "omitted",
            "not implemented",
            "no effect",
            "never",
            "except",
            "without",
        ] {
            assert!(
                !help.contains(word),
                "{word:?} names an exclusion in {help}"
            );
        }
        assert!(!FlagState::implied_by("everything").contains(&"cell_type"));
    }

    /// `--cell_type` is read against the cache's epigenome list by greedy
    /// longest match (`split_cell_types`), since a name may hold a comma:
    /// `A549_(m,_58_y),BJ_(m)` names two epigenomes, and
    /// `common_myeloid_progenitor,_CD34-positive_(27-36_y)` one, where Ensembl
    /// VEP's comma split (Config.pm:761) would ask for `A549_(m` and `_58_y)`.
    /// Trailing commas name nothing more, as Perl's `split` drops the trailing
    /// empty fields, one or several; a name given twice is listed twice.
    #[test]
    fn cell_type_splits_by_known_names() {
        let available: Vec<String> = [
            "A549_(m,_58_y)",
            "BJ_(m)",
            "H1_(m)",
            "H1_(m)_derived_neuronal_progenitor",
            "common_myeloid_progenitor,_CD34-positive_(27-36_y)",
            "astrocyte",
        ]
        .map(str::to_string)
        .to_vec();
        let names = |value: &str| -> Vec<&str> {
            split_cell_types(value, &available)
                .unwrap()
                .into_iter()
                .map(|index| available[index].as_str())
                .collect()
        };
        assert_eq!(names("A549_(m,_58_y),BJ_(m)"), ["A549_(m,_58_y)", "BJ_(m)"]);
        assert_eq!(
            names("common_myeloid_progenitor,_CD34-positive_(27-36_y)"),
            ["common_myeloid_progenitor,_CD34-positive_(27-36_y)"]
        );
        assert_eq!(
            names("H1_(m),H1_(m)_derived_neuronal_progenitor,astrocyte,"),
            ["H1_(m)", "H1_(m)_derived_neuronal_progenitor", "astrocyte"]
        );
        assert_eq!(names("astrocyte,astrocyte"), ["astrocyte", "astrocyte"]);
        assert_eq!(names("astrocyte,,"), ["astrocyte"]);
        assert_eq!(names("A549_(m,_58_y),,,"), ["A549_(m,_58_y)"]);
        assert_eq!(names(""), [] as [&str; 0]);
        assert_eq!(names(",,"), [] as [&str; 0]);
    }

    /// A remainder that starts with no name of the list is refused with Ensembl
    /// VEP's text (AnnotationType/RegFeat.pm:85-92): the piece up to the next
    /// comma, then the list on its own line; a prefix of a name (`astro`), a
    /// name run into the next without a comma (`astrocyteBJ_(m)`) and an empty
    /// piece before a name (`,astrocyte`, `astrocyte,,HUVEC`), which Perl's
    /// `split` keeps where it drops a trailing one, are refused the same way.
    #[test]
    fn cell_type_unknown_name_fails_with_ensembl_veps_text() {
        let available: Vec<String> = ["A549_(m,_58_y)", "BJ_(m)", "HUVEC", "astrocyte"]
            .map(str::to_string)
            .to_vec();
        let refused = |value: &str| split_cell_types(value, &available).unwrap_err().to_string();
        assert_eq!(
            refused("HUVEC,bogus"),
            "ERROR: [cli] Cell type \"bogus\" unavailable; available cell types are:\n\
             A549_(m,_58_y),BJ_(m),HUVEC,astrocyte"
        );
        assert!(refused("bogus,HUVEC").contains("Cell type \"bogus\" unavailable"));
        assert!(refused("astro").contains("Cell type \"astro\" unavailable"));
        assert!(refused("astrocyteBJ_(m)").contains("Cell type \"astrocyteBJ_(m)\" unavailable"));
        assert!(refused(",astrocyte").contains("Cell type \"\" unavailable"));
        assert!(refused("astrocyte,,HUVEC").contains("Cell type \"\" unavailable"));
        assert!(refused("A549_(m,_58_y").contains("Cell type \"A549_(m\" unavailable"));
    }

    /// `--cell_type` switches `--regulatory` on through the option set
    /// (Config.pm:427-432), and its value is kept as given for the cache-time
    /// split; a value Ensembl VEP's comma split reads as an empty list (`""`,
    /// `,`) is not an active flag (`_is_flag_active`, Config.pm:902-905), so it
    /// switches nothing on and names no epigenome.
    #[test]
    fn cell_type_switches_regulatory_on() {
        let c = config(&["vep", "--cell_type", "A549_(m,_58_y),BJ_(m)"]);
        assert!(c.regulatory);
        assert_eq!(c.cell_type.as_deref(), Some("A549_(m,_58_y),BJ_(m)"));
        assert_eq!(FlagState::implied_by("cell_type"), ["regulatory"]);
        for value in ["", ","] {
            let c = config(&["vep", "--cell_type", value]);
            assert!(!c.regulatory, "{value:?}");
            assert!(c.cell_type.is_none(), "{value:?}");
        }
        let c = config(&["vep", "--regulatory"]);
        assert!(c.cell_type.is_none());
    }

    /// `--regulatory_gff` takes Ensembl VEP's grammar (`AnnotationSourceAdaptor.pm:237-262`):
    /// a bare name is the `file`; otherwise comma-separated `key=value` pairs,
    /// refused with its four texts when a pair has no `=` or an empty value, a
    /// key is not `emars`, `file` or `motifs`, or every value is false in Perl.
    /// A trailing comma is no pair, as Perl's `split` drops it; a later pair
    /// overwrites an earlier one with the same key; a value false in Perl
    /// names no file; and the flag's option set switches `--regulatory` on,
    /// except for a value that is itself false.
    #[test]
    fn regulatory_gff_grammar_matches_ensembl_vep() {
        let gff = |value: &str| config(&["vep", "--regulatory_gff", value]);
        let c = gff("reg.gff3.gz");
        assert!(c.regulatory);
        assert_eq!(
            c.regulatory_gff,
            Some(RegulatoryGff {
                file: Some("reg.gff3.gz".into()),
                ..Default::default()
            })
        );
        let c = gff("motifs=m.gff3.gz,file=r.gff3.gz,emars=e.gff.gz,");
        assert_eq!(
            c.regulatory_gff,
            Some(RegulatoryGff {
                file: Some("r.gff3.gz".into()),
                emars: Some("e.gff.gz".into()),
                motifs: Some("m.gff3.gz".into()),
            })
        );
        assert_eq!(
            c.regulatory_gff
                .as_ref()
                .unwrap()
                .sources()
                .collect::<Vec<_>>(),
            [
                ("file", "r.gff3.gz"),
                ("emars", "e.gff.gz"),
                ("motifs", "m.gff3.gz")
            ]
        );
        assert_eq!(
            gff("file=a.gff3.gz,file=b.gff3.gz,motifs=0").regulatory_gff,
            Some(RegulatoryGff {
                file: Some("b.gff3.gz".into()),
                ..Default::default()
            })
        );
        // A bare name that holds a comma is one file name.
        assert_eq!(
            gff("a.gff,b.gff").regulatory_gff.unwrap().file.as_deref(),
            Some("a.gff,b.gff")
        );
        let c = gff("0");
        assert!(c.regulatory_gff.is_none() && !c.regulatory);

        let refused = |value: &str| refusal(&["vep", "--regulatory_gff", value]);
        assert_eq!(
            refused("file=r.gff3.gz,motifs"),
            "ERROR: [cli] Failed to parse --regulatory_gff parameter \"motifs\"; expected key=value"
        );
        assert_eq!(
            refused("file=,motifs=m.gff3.gz"),
            "ERROR: [cli] Failed to parse --regulatory_gff parameter \"file=\"; expected key=value"
        );
        assert_eq!(
            refused("file=r.gff3.gz,,motifs=m.gff3.gz"),
            "ERROR: [cli] Failed to parse --regulatory_gff parameter \"\"; expected key=value"
        );
        assert_eq!(
            refused("matrices=x.json,file=r.gff3.gz,activity=a.tsv"),
            "ERROR: [cli] Unsupported --regulatory_gff key(s): activity, matrices\nSupported keys: emars, file, motifs"
        );
        assert_eq!(
            refused("file=0,emars=0"),
            "ERROR: [cli] --regulatory_gff requires at least one of: emars, file, motifs"
        );
    }

    /// `--extended_promoters` is read by the GFF regulatory source alone
    /// (`RegFeat.pm:170`, a shortcut of that source; no `REQUIRES` entry names
    /// it, `Config.pm:607-612`), so without `--regulatory_gff` the configuration
    /// is the one the same line resolves to without the flag, the flag's own
    /// field and the recorded command line apart.
    #[test]
    fn extended_promoters_without_regulatory_gff_changes_nothing() {
        let with = config(&["vep", "--extended_promoters", "--regulatory"]);
        let without = config(&["vep", "--regulatory"]);
        assert!(with.extended_promoters && !without.extended_promoters);
        assert!(with.regulatory_gff.is_none());
        let settled = Config {
            extended_promoters: false,
            command_line: without.command_line.clone(),
            ..with
        };
        assert_eq!(format!("{settled:?}"), format!("{without:?}"));
        let alone = config(&["vep", "--extended_promoters"]);
        assert!(!alone.regulatory, "the flag switches nothing on");
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

    /// `--parquet_partition` takes `chrom` and nothing else, as a clap value
    /// enumeration, and is refused under the output component with every format
    /// but Parquet, whose configuration carries the value.
    #[test]
    fn parquet_partition_takes_chrom_alone() {
        assert_eq!(
            usage_kind(&["vep", "--parquet_partition", "gene"]),
            clap::error::ErrorKind::InvalidValue
        );
        let message = usage_error(&["vep", "--parquet_partition", "gene"]);
        assert!(
            message.contains("invalid value 'gene' for '--parquet_partition <COLUMN>'")
                && message.contains("[possible values: chrom]"),
            "{message}"
        );
        for argv in [
            &["vep", "--tab", "--parquet_partition", "chrom"][..],
            &["vep", "--parquet_partition", "chrom"][..],
            &["vep", "--vcf", "--parquet_partition", "chrom"][..],
        ] {
            assert_eq!(
                refusal(argv),
                "ERROR: [output] --parquet_partition applies to --output_format parquet",
                "{argv:?}"
            );
        }
        let c = config(&[
            "vep",
            "--output_format",
            "parquet",
            "--parquet_partition",
            "chrom",
            "-o",
            "x.parquet",
        ]);
        assert_eq!(c.parquet_partition.as_deref(), Some("chrom"));
        assert_eq!(
            config(&["vep", "--output_format", "parquet", "-o", "x.parquet"]).parquet_partition,
            None
        );
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

    /// `--all_refseq` (`Config.pm:98`) is accepted beside `--refseq` and beside
    /// `--merged`, the caches whose RefSeq clause it switches off
    /// (`AnnotationType/Transcript.pm:237`), and alone, where Ensembl VEP accepts
    /// it too and the clause never runs; `%INCOMPATIBLE` (`Config.pm:615-631`)
    /// names no pair with it.
    #[test]
    fn all_refseq_is_accepted_beside_refseq_and_merged() {
        use crate::transcript_index::TranscriptSet;

        let refseq = config(&["vep", "--refseq", "--all_refseq"]);
        assert!(refseq.refseq && refseq.all_refseq);
        let merged = config(&["vep", "--merged", "--all_refseq"]);
        assert!(merged.merged && merged.all_refseq);
        let alone = config(&["vep", "--all_refseq"]);
        assert!(alone.all_refseq && !alone.refseq && !alone.merged);
        assert_eq!(
            TranscriptSet::from_config(&alone),
            TranscriptSet {
                all_refseq: true,
                ..TranscriptSet::default()
            }
        );
        assert!(TranscriptSet::from_config(&alone).accepts_model(
            "LOC100419737",
            "21",
            None,
            false,
            false
        ));
        assert!(!config(&["vep"]).all_refseq);
        assert!(FlagState::implied_by("all_refseq").is_empty());
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

    /// A `--convert` value Perl reads as false, `0` or the empty value, leaves
    /// the flag inactive (Config.pm:902-906), so the deprecation check passes it
    /// (:977) and the run proceeds; any other value is refused, `00` included.
    #[test]
    fn convert_with_a_false_value_proceeds() {
        for argv in [
            &["vep", "--convert", "0"][..],
            &["vep", "--convert", ""][..],
            &["vep", "--convert=0"][..],
            &["vep", "--convert="][..],
        ] {
            assert!(Config::from_argv(argv).is_ok(), "{argv:?}");
        }
        assert_eq!(
            error_table_entry("convert").unwrap().help(),
            "Refused with an error unless its value is 0 or empty: --convert has been deprecated"
        );
        for value in ["vcf", "00", "1"] {
            assert_eq!(
                refusal(&["vep", "--convert", value]),
                "ERROR: [cli] --convert has been deprecated",
                "{value}"
            );
        }
    }

    /// `custom` is read from the state after the option sets, which the `gff`,
    /// `gtf`, `bigwig`, `phyloP` and `phastCons` sets append to (Config.pm:531-564);
    /// `gff` and `gtf` are refused before the sets run, so `bigwig` stands in.
    #[test]
    fn custom_follows_the_state_the_option_sets_append_to() {
        let mut parsed = parse_command_line(["vep", "--custom", "x.bed,X,bed"]).unwrap();
        parsed.state.set("bigwig", "a.bw");
        let c = Config::from_parsed(parsed).unwrap();
        assert_eq!(c.custom, ["x.bed,X,bed", "a.bw,,bigwig,exact"]);
        assert_eq!(
            config(&["vep", "--custom", "y.bed,Y,bed"]).custom,
            ["y.bed,Y,bed"]
        );
        assert!(config(&["vep"]).custom.is_empty());
    }

    /// `--custom_multi_allelic` needs a `--custom` source (`%REQUIRES`,
    /// Config.pm:611; the refusal text Config.pm:972), and is read with one.
    #[test]
    fn custom_multi_allelic_requires_custom() {
        assert_eq!(
            refusal(&["vep", "--custom_multi_allelic"]),
            "ERROR: [cli] You must set --custom to use --custom_multi_allelic"
        );
        let c = config(&[
            "vep",
            "--custom",
            "x.vcf.gz,X,vcf",
            "--custom_multi_allelic",
        ]);
        assert!(c.custom_multi_allelic);
        assert!(!config(&["vep", "--custom", "x.vcf.gz,X,vcf"]).custom_multi_allelic);
    }

    /// `--custom_suppress_filter` is a switch Ensembl VEP declares without a
    /// requirement (Config.pm:268).
    #[test]
    fn custom_suppress_filter_is_accepted() {
        assert!(config(&["vep", "--custom_suppress_filter"]).custom_suppress_filter);
        assert!(!config(&["vep"]).custom_suppress_filter);
    }

    /// `--bigwig <file>` appends the positional source `<file>,,bigwig,exact`
    /// (`Config.pm:545-550`), whose short name is the file's last path segment
    /// (`File.pm:154-157`), after the command line's own `--custom` sources.
    #[test]
    fn bigwig_shortcut_appends_a_positional_custom_source() {
        assert_eq!(
            config(&["vep", "--bigwig", "x.bw"]).custom,
            ["x.bw,,bigwig,exact"]
        );
        let spec = vep_builtins::custom::parse_custom("x.bw,,bigwig,exact").unwrap();
        assert_eq!(spec.short_name, "x.bw");
        assert_eq!(spec.match_type, vep_builtins::custom::MatchType::Exact);
        assert_eq!(
            config(&["vep", "--custom", "a.bed,A,bed", "--bigwig", "/d/x.bw"]).custom,
            ["a.bed,A,bed", "/d/x.bw,,bigwig,exact"]
        );
    }

    /// Each `--phyloP <N>` appends one bigWig source under the local tree's layout
    /// (`Config.pm:552-557`, one per instance, `:829-835`), its short name spelt
    /// `phyloP<N>way`: Ensembl VEP's set spells it `phlyoP<N>way`
    /// (`Config.pm:555`), a misspelling whose column vep-rs does not print.
    #[test]
    fn phylop_expands_one_bigwig_source_per_instance_with_the_corrected_short_name() {
        let c = config(&[
            "vep",
            "--phyloP",
            "100",
            "--phyloP",
            "20",
            "--ucsc_assembly",
            "hg38",
            "--ucsc_data_root",
            "/m/",
        ]);
        assert_eq!(
            c.custom,
            [
                "/m/hg38/phyloP100way/hg38.phyloP100way.bw,phyloP100way,bigwig,exact",
                "/m/hg38/phyloP20way/hg38.phyloP20way.bw,phyloP20way,bigwig,exact",
            ]
        );
        assert!(!c.custom.iter().any(|source| source.contains("phlyoP")));
    }

    /// `--phyloP` without `--ucsc_assembly` fails with Ensembl VEP's text
    /// (`Config.pm:609, 972`), before the root is looked at.
    #[test]
    fn phylop_without_ucsc_assembly_fails_with_vep_text() {
        assert_eq!(
            refusal(&["vep", "--phyloP", "100"]),
            "ERROR: [cli] You must set --ucsc_assembly to use --phyloP"
        );
        assert_eq!(
            refusal(&["vep", "--phastCons", "100", "--ucsc_data_root", "/m/"]),
            "ERROR: [cli] You must set --ucsc_assembly to use --phastCons"
        );
    }

    /// Each `--phastCons <N>` appends one bigWig source under the local tree's
    /// layout (`Config.pm:559-564`), its short name `phastCons<N>way` as
    /// Ensembl VEP spells it.
    #[test]
    fn phastcons_expands_one_bigwig_source_per_instance() {
        let c = config(&[
            "vep",
            "--phastCons",
            "100",
            "--phastCons",
            "30",
            "--ucsc_assembly",
            "hg19",
            "--ucsc_data_root",
            "/m/",
        ]);
        assert_eq!(
            c.custom,
            [
                "/m/hg19/phastCons100way/hg19.phastCons100way.bw,phastCons100way,bigwig,exact",
                "/m/hg19/phastCons30way/hg19.phastCons30way.bw,phastCons30way,bigwig,exact",
            ]
        );
    }

    /// `--ucsc_assembly` stands in both path segments of a shortcut's file, the
    /// directory and the file name (`Config.pm:555, 562`); any name is accepted
    /// (`Config.pm:596-604` lists no values for it).
    #[test]
    fn ucsc_assembly_is_substituted_into_both_path_segments() {
        let c = config(&[
            "vep",
            "--phyloP",
            "7",
            "--phastCons",
            "7",
            "--ucsc_assembly",
            "mm39",
            "--ucsc_data_root",
            "/mirror/",
        ]);
        assert_eq!(
            c.custom,
            [
                "/mirror/mm39/phyloP7way/mm39.phyloP7way.bw,phyloP7way,bigwig,exact",
                "/mirror/mm39/phastCons7way/mm39.phastCons7way.bw,phastCons7way,bigwig,exact",
            ]
        );
    }

    /// The default root is Ensembl VEP's remote one (`Config.pm:314`), which a
    /// shortcut refuses with the local tree's layout and the `--custom` form of the
    /// same source; without a shortcut the default is inert, and an explicit
    /// remote root is refused the same way.
    #[test]
    fn ucsc_data_root_default_is_refused_with_a_local_mirror_message() {
        let message = "ERROR: [custom] --ucsc_data_root must name a local directory laid out as \
                       <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw; the same \
                       source is --custom file=<local>.bw,format=bigwig,type=exact";
        assert_eq!(
            refusal(&["vep", "--phyloP", "100", "--ucsc_assembly", "hg38"]),
            message
        );
        assert_eq!(
            refusal(&[
                "vep",
                "--phastCons",
                "100",
                "--ucsc_assembly",
                "hg38",
                "--ucsc_data_root",
                "ftp://mirror.example.org/goldenpath/"
            ]),
            message
        );
        assert!(config(&["vep"]).custom.is_empty());
        assert!(config(&["vep", "--ucsc_data_root", "http://example.org/"])
            .custom
            .is_empty());
    }

    /// A local `--ucsc_data_root` is joined to the assembly as given
    /// (`Config.pm:555`: `%ucsc_data_root%%ucsc_assembly%`, no separator), so
    /// the directory is named with its trailing slash.
    #[test]
    fn ucsc_data_root_local_directory_builds_the_bigwig_path() {
        let c = config(&[
            "vep",
            "--phyloP",
            "100",
            "--ucsc_assembly",
            "hg38",
            "--ucsc_data_root",
            "/data/ucsc/",
        ]);
        assert_eq!(
            c.custom,
            ["/data/ucsc/hg38/phyloP100way/hg38.phyloP100way.bw,phyloP100way,bigwig,exact"]
        );
        assert_eq!(
            config(&[
                "vep",
                "--phyloP",
                "100",
                "--ucsc_assembly",
                "hg38",
                "--ucsc_data_root",
                "mirror"
            ])
            .custom,
            ["mirrorhg38/phyloP100way/hg38.phyloP100way.bw,phyloP100way,bigwig,exact"]
        );
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

    /// Ensembl VEP 116.2 reads `-i` as standard input when its value is `STDIN`
    /// in any letter case (`uc($file) eq 'STDIN'`, Parser.pm:235) and keeps the
    /// parameter as given, so `-i stdin` prints `Input file<TAB>stdin` in the
    /// statistics report (Stats.pm:652) and `stdin` in the skipped-variants
    /// file's heading (BaseRunner.pm:371-372); the default is `STDIN`
    /// (Config.pm:301).
    #[test]
    fn input_file_keeps_the_stdin_spelling_given() {
        for spelling in ["stdin", "Stdin", "STDIN"] {
            let config = config(&["vep", "-i", spelling]);
            assert_eq!(config.input_file, spelling);
            assert!(config.reads_stdin(), "{spelling}");
        }
        assert_eq!(config(&["vep"]).input_file, "STDIN");
        assert!(config(&["vep"]).reads_stdin());
        let named = config(&["vep", "-i", "in.vcf"]);
        assert_eq!(named.input_file, "in.vcf");
        assert!(!named.reads_stdin());
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
            "parquet",
            "cache",
            "fasta",
            "filter",
            "check_svs",
            "custom",
            "config",
            "run",
        ] {
            assert!(names.contains(component), "{component}");
        }
    }

    /// `--no_check_variants_order` is a flag of its own, outside the table of
    /// flags accepted without effect, and reaches the configuration as given.
    #[test]
    fn no_check_variants_order_reaches_the_config() {
        assert!(error_table_entry("no_check_variants_order").is_none());
        assert!(!config(&["vep"]).no_check_variants_order);
        assert!(config(&["vep", "--no_check_variants_order"]).no_check_variants_order);
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

    /// `shift_3prime=i` and `shift_genomic=i` (Config.pm:231-232) are read for
    /// truth alone (AnnotationType/Transcript.pm:162, OutputFactory.pm:1349), on
    /// the number Getopt::Long stores: any integer but 0 is on, and a spelling
    /// of zero is off and not active for the pair `--shift_hgvs` refuses
    /// (Config.pm:630).
    #[test]
    fn shift_3prime_and_shift_genomic_are_on_for_any_integer_but_zero() {
        for value in [
            "1",
            "2",
            "10",
            "-1",
            "+1",
            "01",
            "1_0",
            "9999999999999999999999",
        ] {
            let c = config(&["vep", "--shift_3prime", value, "--shift_genomic", value]);
            assert!(c.shift_3prime && c.shift_genomic, "{value}");
            for flag in ["--shift_3prime", "--shift_genomic"] {
                assert_eq!(
                    refusal(&["vep", "--shift_hgvs", "1", flag, value]),
                    format!("ERROR: [cli] Can't use --shift_hgvs and {flag} together"),
                    "{flag} {value}"
                );
            }
        }
        for value in ["0", "00", "-0", "+0", "0_0"] {
            let c = config(&["vep", "--shift_3prime", value, "--shift_genomic", value]);
            assert!(!c.shift_3prime && !c.shift_genomic, "{value}");
            let with_shift_hgvs = [
                "vep",
                "--shift_hgvs",
                "1",
                "--shift_3prime",
                value,
                "--shift_genomic",
                value,
            ];
            assert!(Config::from_argv(with_shift_hgvs).is_ok(), "{value}");
        }
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

    /// The plugins directory value is `--dir_plugins`, else `<dir>/Plugins`
    /// (Runner.pm:907, whose `||` passes an empty string and `0` over as
    /// Perl-false), a leading `~` expanded as `glob` expands it (`:908`); a value
    /// with no metacharacter is its one match whether or not it exists.
    #[test]
    fn dir_plugins_defaults_to_dir_slash_plugins() {
        let home = std::env::var("HOME").unwrap();
        let matches = |argv: &[&str]| -> Vec<String> {
            let dirs = config(argv).dir_plugins;
            let mut all = dirs.not_found;
            all.extend(dirs.dirs);
            all
        };
        assert_eq!(matches(&["vep"]), [format!("{home}/.vep/Plugins")]);
        assert_eq!(
            matches(&["vep", "--dir", "/data/vep"]),
            ["/data/vep/Plugins"]
        );
        assert_eq!(
            matches(&["vep", "--dir", "/data/vep", "--dir_plugins", "/p"]),
            ["/p"]
        );
        for value in ["", "0"] {
            assert_eq!(
                matches(&["vep", "--dir", "/data/vep", "--dir_plugins", value]),
                ["/data/vep/Plugins"],
                "{value:?}"
            );
        }
        assert_eq!(
            matches(&["vep", "--dir_plugins", "~/plugins"]),
            [format!("{home}/plugins")]
        );
        assert_eq!(
            config(&["vep", "--dir_plugins", "~/plugins"])
                .dir_plugins
                .pattern,
            "~/plugins"
        );
        assert!(!config(&["vep"]).safe);
        assert!(config(&["vep", "--safe"]).safe);
    }

    /// `--dir_plugins` is expanded as `glob` expands it (Runner.pm:908): over a
    /// tree of `p1/`, `p2/` and `f.txt`, `<tmp>/*` matches all three; the two
    /// directories are searched last match first, the order `unshift @INC` over
    /// the matches leaves (`:921`), and the file is the one match to warn about
    /// (`:911-916`). A pattern matching nothing names no directory and the
    /// `Failed to load plugin` line falls back to the value as given; a value
    /// naming one missing directory is reported as that directory.
    #[test]
    fn dir_plugins_glob_expands_to_every_directory_match() {
        let tmp = tempfile::tempdir().unwrap();
        for name in ["p1", "p2"] {
            std::fs::create_dir(tmp.path().join(name)).unwrap();
        }
        std::fs::write(tmp.path().join("f.txt"), "").unwrap();
        let root = tmp.path().to_str().unwrap();
        let dirs = config(&["vep", "--dir_plugins", &format!("{root}/*")]).dir_plugins;
        assert_eq!(dirs.dirs, [format!("{root}/p2"), format!("{root}/p1")]);
        assert_eq!(dirs.not_found, [format!("{root}/f.txt")]);
        assert_eq!(dirs.pattern, format!("{root}/*"));
        assert_eq!(dirs.searched(), format!("{root}/p2, {root}/p1"));

        let none = config(&["vep", "--dir_plugins", &format!("{root}/zz*")]).dir_plugins;
        assert_eq!(
            none,
            PluginDirs {
                pattern: format!("{root}/zz*"),
                ..PluginDirs::default()
            }
        );
        assert_eq!(none.searched(), format!("{root}/zz*"));

        let missing = PluginDirs::expand("/no/such/dir");
        assert_eq!(missing.not_found, ["/no/such/dir"]);
        assert!(missing.dirs.is_empty());
        assert_eq!(missing.searched(), "/no/such/dir");
    }

    const DYLIB_EXT: &str = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };

    /// The profile's `examples/` directory when the `echo_plugin` and
    /// `echo_all_rows_plugin` example dylibs have been built there, which a
    /// `cargo test` selection including `vep-plugin` does.
    fn example_plugins_dir() -> Option<std::path::PathBuf> {
        let exe = std::env::current_exe().unwrap();
        let examples = exe.parent().unwrap().parent().unwrap().join("examples");
        let built = ["echo_plugin", "echo_all_rows_plugin"]
            .iter()
            .all(|name| examples.join(format!("lib{name}.{DYLIB_EXT}")).is_file());
        if built {
            return Some(examples);
        }
        assert!(
            std::env::var_os("CI").is_none(),
            "example dylibs not found under {}",
            examples.display()
        );
        eprintln!(
            "Skipping: example dylibs not found under {}",
            examples.display()
        );
        None
    }

    /// Two plugins directories each holding a `libecho_plugin` library: the one
    /// in the match that sorts last loads, since `unshift @INC` over the matches
    /// in order (Runner.pm:921) searches that directory first, and a brace group
    /// orders its members as written. The two example libraries stand in for two
    /// versions of one plugin: `echo_plugin` reports version `0.2.0`,
    /// `echo_all_rows_plugin` `0.1.0`.
    #[test]
    fn dir_plugins_search_order_is_last_match_first() {
        let Some(examples) = example_plugins_dir() else {
            return;
        };
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().to_str().unwrap();
        for (dir, example) in [("p1", "echo_all_rows_plugin"), ("p2", "echo_plugin")] {
            std::fs::create_dir(tmp.path().join(dir)).unwrap();
            std::fs::copy(
                examples.join(format!("lib{example}.{DYLIB_EXT}")),
                tmp.path()
                    .join(dir)
                    .join(format!("libecho_plugin.{DYLIB_EXT}")),
            )
            .unwrap();
        }
        let version_under = |pattern: &str| -> String {
            let dirs = PluginDirs::expand(&format!("{root}/{pattern}"));
            let loaded = crate::annotator::load_plugins(
                &["echo_plugin".to_string()],
                &dirs,
                true,
                None,
                None,
                vep_builtins::traits::RunContext::default(),
            )
            .unwrap();
            loaded.dylib.plugins()[0].version().to_string()
        };
        assert_eq!(version_under("*"), "0.2.0");
        assert_eq!(version_under("{p2,p1}"), "0.1.0");
        assert_eq!(version_under("{p1,p2}"), "0.2.0");
        assert_eq!(version_under("p1"), "0.1.0");
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

    /// `Runner.pm:787-790`: the HTML report is on unless `--stats_text` alone
    /// names a report; `--stats_html` keeps it on beside `--stats_text`, and
    /// `--no_stats` is read separately (Runner.pm:144 tests it first).
    #[test]
    fn stats_html_is_on_unless_stats_text_alone_names_a_report() {
        let default = config(&["vep"]);
        assert!(default.stats_html && !default.stats_text && !default.no_stats);
        assert!(config(&["vep", "--stats_html"]).stats_html);
        let text = config(&["vep", "--stats_text"]);
        assert!(text.stats_text && !text.stats_html);
        let both = config(&["vep", "--stats_text", "--stats_html"]);
        assert!(both.stats_text && both.stats_html);
        let silent = config(&["vep", "--no_stats"]);
        assert!(silent.no_stats && silent.stats_html);
        assert!(config(&["vep", "--stats_file", "s.html"]).stats_html);
    }
}
