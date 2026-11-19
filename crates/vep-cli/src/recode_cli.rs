// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep recode`: the command line of the companion tool that translates a
//! variant between its notations, as Ensembl VEP's `variant_recoder` does
//! (`variant_recoder` at release/116.2, whose line numbers the citations below
//! name). The translation itself is [`crate::recode`].
//!
//! The options are `variant_recoder`'s (36-95) less the ones that name an
//! Ensembl database connection or a RefSeq alignment file, which are refused by
//! name: the tool reads the cache alone. `--refseq` is refused with the
//! recoder's own text, since its forced `merged` (`VariantRecoder.pm:98`) meets
//! `Config.pm:621` on every run. The options are read as Getopt::Long reads
//! them, after one dash or two, in any case, abbreviated to any prefix that
//! names one option ([`crate::argv::spell_as_declared`]); a repeated option
//! keeps its last value, and every `--transcript_filter` given is kept (66).
//! Behind the command line, the `--config` file, `<dir>/vep.ini` and the
//! `VEP_*` variables are read as `Config->new` reads them for the recoder,
//! each lower source yielding to the ones above it (`Config.pm:701-708`,
//! [`crate::argv::with_config_file`]), the options the recoder presets or
//! reads before them left to the command line ([`PRESET_OPTIONS`]).
//! Every error is one line `ERROR: [<component>] <message>` on standard error
//! and exit 1, the component naming what refused the run (`recode`, `config`,
//! `input`, `cache`, `fasta` or `filter`); a command-line syntax error is
//! clap's, exit 2; no argument at all prints the usage and exits 0 (97).

use std::ffi::OsString;
use std::io::{self, Write};

use clap::error::ErrorKind as ClapErrorKind;
use clap::{ArgAction, CommandFactory, Parser};

use crate::argv;
use crate::config::CliError;
use crate::recode;

const COMPONENT: &str = "recode";

/// The options a source behind the command line never sets, as `Config->new`
/// finds them for `variant_recoder`: `fields`, `buffer_size`, `failed`,
/// `merged`, `lrg` and `ambiguous_hgvs` are preset before it runs
/// (`VariantRecoder.pm:96-111`, `:152`), so its `||=` (`Config.pm:1101`,
/// `:1158`) leaves them; `vcf_string`, `var_synonyms`, `mane_select` and
/// `ga4gh_vrs` add their fields at `:131-150`, before the sources are read, so
/// a source's value adds none; `pretty` and `vcf_string` are no options of
/// `@VEP_PARAMS`, which `is_valid_param` (`Config.pm:1046`) drops. Each is set
/// on the command line alone.
pub const PRESET_OPTIONS: &[&str] = &[
    "fields",
    "buffer_size",
    "failed",
    "merged",
    "lrg",
    "ambiguous_hgvs",
    "vcf_string",
    "var_synonyms",
    "mane_select",
    "ga4gh_vrs",
    "pretty",
];

fn error(message: impl Into<String>) -> CliError {
    CliError::new(COMPONENT, message)
}

/// The command line of `vep recode`, `variant_recoder`'s options (36-95),
/// `--config` (41) and `--debug` (94) among them. A repeated option overrides
/// itself, as Getopt::Long assigns each occurrence of a scalar option over the
/// last (`Getopt/Long.pm:677`); `--transcript_filter` appends, as its array
/// destination does (66).
#[derive(Debug, Parser)]
#[command(
    name = "vep recode",
    about = "Translate a variant between its identifiers and notations, as Ensembl VEP's variant_recoder does",
    long_about = "Translate a variant between its identifiers and notations, as Ensembl VEP's variant_recoder does.\n\nThe input is an identifier, an HGVS, SPDI, VCF, region or Ensembl-format record, one per line, read from --input_file or given as --input_data. Each input line becomes one object of a JSON array, keyed by the input's alleles, whose values hold the fields named as arrays: id, hgvsg, hgvsc, hgvsp and spdi by default, plus vcf_string, var_synonyms, mane_select and ga4gh_vrs under their flags. Every field comes from the cache alone.",
    disable_version_flag = true,
    arg_required_else_help = false,
    args_override_self = true
)]
pub struct RecodeArgs {
    /// A configuration file, one option per line as <name> <value> without the dashes (cache_version 116, pick 1); an option the command line also gives keeps the command line's value
    #[arg(long = "config", value_name = "FILE", help_heading = "Input")]
    pub config: Option<String>,

    /// The input file, one variant per line; a gzip or BGZF file is decompressed
    #[arg(
        long = "input_file",
        short = 'i',
        value_name = "FILE",
        help_heading = "Input"
    )]
    pub input_file: Option<String>,

    /// The input as a string, read in place of --input_file; one variant per line
    #[arg(
        long = "input_data",
        visible_alias = "id",
        value_name = "TEXT",
        help_heading = "Input"
    )]
    pub input_data: Option<String>,

    /// The input format: ensembl, vcf, hgvs, id, spdi or region, detected from the first line when absent
    #[arg(long = "format", value_name = "FORMAT", help_heading = "Input")]
    pub format: Option<String>,

    /// The column delimiter of an Ensembl-format or VCF input
    #[arg(long = "delimiter", value_name = "TEXT", help_heading = "Input")]
    pub delimiter: Option<String>,

    /// Indent the JSON output, three spaces per level, as Ensembl VEP prints it
    #[arg(long = "pretty", help_heading = "Output")]
    pub pretty: bool,

    /// The fields of each allele's object, comma-separated; default id,hgvsg,hgvsc,hgvsp,spdi. A field named here switches on what computes it: hgvsg, hgvsc, hgvsp and spdi need the FASTA, id the co-located variants. Any key of the annotation's JSON output may be named
    #[arg(long = "fields", value_name = "LIST", help_heading = "Output")]
    pub fields: Option<String>,

    /// Add vcf_string: each allele's VCF form as chromosome-position-reference-alternate
    #[arg(long = "vcf_string", help_heading = "Output")]
    pub vcf_string: bool,

    /// Add var_synonyms: the known synonyms of the variant, one "source: id,id" string per source, from the cache's co-located variants
    #[arg(long = "var_synonyms", help_heading = "Output")]
    pub var_synonyms: bool,

    /// Add mane_select: the hgvsg, hgvsc and hgvsp of the MANE Select transcript, which switches those three fields on; human only
    #[arg(long = "mane_select", help_heading = "Output")]
    pub mane_select: bool,

    /// Add ga4gh_vrs: a GA4GH VRS Allele object for each spdi, which it switches on; an empty list for a chromosome whose synonyms name no NC_ accession
    #[arg(long = "ga4gh_vrs", help_heading = "Output")]
    pub ga4gh_vrs: bool,

    /// Also write the warnings, which the output carries, to this file
    #[arg(long = "warning_file", value_name = "FILE", help_heading = "Output")]
    pub warning_file: Option<String>,

    /// Accepted as variant_recoder accepts it: there it prints the plugin directories and modules loaded, of which vep recode has none, so the output is unchanged
    #[arg(long = "debug", help_heading = "Output")]
    pub debug: bool,

    /// The cache directory itself, the one holding info.json and transcripts/, read in place of the --dir_cache layout
    #[arg(
        long = "json_cache",
        value_name = "DIR",
        help_heading = "Cache, species and assembly"
    )]
    pub json_cache: Option<String>,

    /// The cache root: the cache read is <dir_cache>/<species>/<cache_version>_<assembly>, the species directory suffixed _merged when the merged set is read, as Ensembl VEP lays it out
    #[arg(
        long = "dir_cache",
        value_name = "DIR",
        help_heading = "Cache, species and assembly"
    )]
    pub dir_cache: Option<String>,

    /// The cache version, the prefix of the cache directory under the species directory
    #[arg(
        long = "cache_version",
        value_name = "N",
        help_heading = "Cache, species and assembly"
    )]
    pub cache_version: Option<u32>,

    /// The reference FASTA, needed for the HGVS and SPDI fields; a FASTA found in the cache directory is read without the flag
    #[arg(
        long = "fasta",
        value_name = "FILE",
        help_heading = "Cache, species and assembly"
    )]
    pub fasta: Option<String>,

    /// A chromosome synonyms file, read in place of the cache's own
    #[arg(
        long = "synonyms",
        value_name = "FILE",
        help_heading = "Cache, species and assembly"
    )]
    pub synonyms: Option<String>,

    /// The species, as the cache directory names it; default homo_sapiens
    #[arg(
        long = "species",
        short = 's',
        value_name = "SPECIES",
        help_heading = "Cache, species and assembly"
    )]
    pub species: Option<String>,

    /// The assembly, checked against the cache's
    #[arg(
        long = "assembly",
        short = 'a',
        value_name = "NAME",
        help_heading = "Cache, species and assembly"
    )]
    pub assembly: Option<String>,

    /// Select the GRCh37 assembly
    #[arg(long = "grch37", help_heading = "Cache, species and assembly")]
    pub grch37: bool,

    /// Refused: the recoder reads the merged transcript set on every run, so the flag meets "Can't use --refseq and --merged together"
    #[arg(long = "refseq", help_heading = "Transcript selection")]
    pub refseq: bool,

    /// Read the merged transcript set, Ensembl and RefSeq models together, the species directory suffixed _merged; a cache that holds the set is read without the flag, as variant_recoder presets it
    #[arg(long = "merged", help_heading = "Transcript selection")]
    pub merged: bool,

    /// Keep every model of a merged cache, not the RefSeq accessions alone
    #[arg(long = "all_refseq", help_heading = "Transcript selection")]
    pub all_refseq: bool,

    /// Use the GENCODE basic transcript set alone
    #[arg(long = "gencode_basic", help_heading = "Transcript selection")]
    pub gencode_basic: bool,

    /// A transcript filter, as vep's --transcript_filter takes it; repeatable
    #[arg(
        long = "transcript_filter",
        value_name = "FILTER",
        action = ArgAction::Append,
        allow_hyphen_values = true,
        help_heading = "Transcript selection"
    )]
    pub transcript_filter: Vec<String>,

    /// Drop RefSeq's predicted transcript models, the XM_ and XR_ accessions
    #[arg(long = "exclude_predicted", help_heading = "Transcript selection")]
    pub exclude_predicted: bool,

    /// Convert each input allele to its minimal representation before annotation
    #[arg(long = "minimal", help_heading = "Transcript selection")]
    pub minimal: bool,

    /// Keep one consequence per variant
    #[arg(long = "pick", help_heading = "Transcript selection")]
    pub pick: bool,

    /// Keep one consequence per variant allele
    #[arg(long = "pick_allele", help_heading = "Transcript selection")]
    pub pick_allele: bool,

    /// Keep one consequence per gene
    #[arg(long = "per_gene", help_heading = "Transcript selection")]
    pub per_gene: bool,

    /// Keep one consequence per allele per gene
    #[arg(long = "pick_allele_gene", help_heading = "Transcript selection")]
    pub pick_allele_gene: bool,

    /// The criteria the pick flags compare, in order, comma-separated
    #[arg(
        long = "pick_order",
        value_name = "LIST",
        help_heading = "Transcript selection"
    )]
    pub pick_order: Option<String>,

    /// The number of variants annotated per batch; default 1
    #[arg(
        long = "buffer_size",
        value_name = "N",
        help_heading = "Transcript selection"
    )]
    pub buffer_size: Option<usize>,

    /// The highest failed value a known variant may carry and still be reported; default 1
    #[arg(
        long = "failed",
        value_name = "N",
        help_heading = "Transcript selection"
    )]
    pub failed: Option<i64>,

    /// Describe an insertion or deletion at its most 3' position in HGVS (1, the default) or where the input placed it (0)
    #[arg(
        long = "shift_hgvs",
        value_name = "0|1",
        help_heading = "Transcript selection"
    )]
    pub shift_hgvs: Option<String>,

    /// Read each allele's reference from the transcript instead of the input
    #[arg(long = "use_transcript_ref", help_heading = "Transcript selection")]
    pub use_transcript_ref: bool,

    /// Keep the input's reference allele on a cache that would otherwise switch --use_transcript_ref on
    #[arg(long = "use_given_ref", help_heading = "Transcript selection")]
    pub use_given_ref: bool,

    /// Refused: variant_recoder reads an Ensembl database; vep recode reads a local cache alone
    #[arg(long = "registry", value_name = "FILE", help_heading = "Database")]
    pub registry: Option<String>,

    /// Refused: the database host
    #[arg(long = "host", value_name = "HOST", help_heading = "Database")]
    pub host: Option<String>,

    /// Refused: the database port
    #[arg(long = "port", value_name = "PORT", help_heading = "Database")]
    pub port: Option<String>,

    /// Refused: the database user
    #[arg(
        long = "user",
        short = 'u',
        value_name = "USER",
        help_heading = "Database"
    )]
    pub user: Option<String>,

    /// Refused: the database password
    #[arg(
        long = "password",
        visible_alias = "pass",
        value_name = "PASS",
        help_heading = "Database"
    )]
    pub password: Option<String>,

    /// Refused: the Ensembl database version; --cache_version selects a cache
    #[arg(long = "db_version", value_name = "N", help_heading = "Database")]
    pub db_version: Option<String>,

    /// Refused: the Ensembl Genomes database host
    #[arg(long = "genomes", help_heading = "Database")]
    pub genomes: bool,

    /// Refused: a multi-species Ensembl database
    #[arg(
        long = "is_multispecies",
        value_name = "0|1",
        help_heading = "Database"
    )]
    pub is_multispecies: Option<String>,

    /// Refused: LRG mapping needs the Ensembl core database
    #[arg(long = "lrg", help_heading = "Database")]
    pub lrg: bool,

    /// Refused: resolving a protein notation to every matching variant needs the Ensembl database
    #[arg(long = "ambiguous_hgvs", help_heading = "Database")]
    pub ambiguous_hgvs: bool,

    /// Refused: an alignment file that would edit the cache's transcript models
    #[arg(long = "bam", value_name = "FILE", help_heading = "Database")]
    pub bam: Option<String>,
}

impl RecodeArgs {
    /// The flags the tool refuses, by name, with the reason: the ones that name
    /// an Ensembl database connection (`variant_recoder:48-63`), the database-only
    /// features and the RefSeq alignment file (`:89`).
    pub fn refused_flag(&self) -> Option<CliError> {
        let database = "names the Ensembl database connection variant_recoder reads; vep recode \
                        reads a local cache alone, named with --json_cache or --dir_cache";
        let refusals: [(&str, bool, &str); 11] = [
            ("registry", self.registry.is_some(), database),
            ("host", self.host.is_some(), database),
            ("port", self.port.is_some(), database),
            ("user", self.user.is_some(), database),
            ("password", self.password.is_some(), database),
            (
                "db_version",
                self.db_version.is_some(),
                "selects an Ensembl database version; --cache_version selects a cache",
            ),
            (
                "genomes",
                self.genomes,
                "selects the Ensembl Genomes database host; vep recode reads a local cache alone",
            ),
            (
                "is_multispecies",
                self.is_multispecies.is_some(),
                "describes a multi-species Ensembl database; vep recode reads a local cache alone",
            ),
            (
                "lrg",
                self.lrg,
                "maps a variant to its LRG transcripts through the Ensembl core database, which vep recode does not read",
            ),
            (
                "ambiguous_hgvs",
                self.ambiguous_hgvs,
                "resolves a protein notation to every matching variant through the Ensembl database; a protein notation resolves to one variant here",
            ),
            (
                "bam",
                self.bam.is_some(),
                "names an alignment file that would edit the cache's transcript models; a cache built against one carries the edits itself",
            ),
        ];
        refusals
            .iter()
            .find(|(_, given, _)| *given)
            .map(|(flag, _, reason)| error(format!("--{flag} is not supported: it {reason}")))
    }
}

/// Runs `vep recode` over `argv`, the tokens after the `recode` word, and
/// returns the exit status: 0 on success, 1 after an `ERROR: [<component>]`
/// line, 2 on a command-line syntax error.
pub fn run(argv: Vec<OsString>) -> i32 {
    let mut command = RecodeArgs::command();
    command.build();
    if argv.is_empty() {
        // variant_recoder prints its usage and exits 0 with no argument (97).
        return if command.print_long_help().is_ok() {
            0
        } else {
            1
        };
    }
    let spelled = argv::spell_as_declared(&argv, &command);
    // `Config->new` reads the `--config` file, `vep.ini` and the `VEP_*`
    // variables behind the command line (Config.pm:701-708), so a flag they
    // name reaches the parser as one the command line gave.
    let spelled =
        match argv::with_config_file(&command, &spelled, std::env::vars_os(), PRESET_OPTIONS) {
            Ok(spelled) => spelled,
            Err(error) => {
                eprintln!("{error}");
                return 1;
            }
        };
    let args = match RecodeArgs::try_parse_from(
        std::iter::once(OsString::from("vep recode")).chain(spelled),
    ) {
        Ok(args) => args,
        Err(error) if error.kind() == ClapErrorKind::DisplayHelp => {
            let _ = error.print();
            return 0;
        }
        Err(error) => error.exit(),
    };
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    match recode::run(args, &mut out).and_then(|()| out.flush().map_err(|e| error(e.to_string()))) {
        Ok(()) => 0,
        Err(error) => {
            eprintln!("{error}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(argv: &[&str]) -> RecodeArgs {
        let mut command = RecodeArgs::command();
        command.build();
        let argv: Vec<OsString> = argv.iter().map(OsString::from).collect();
        let spelled = argv::spell_as_declared(&argv, &command);
        RecodeArgs::try_parse_from(std::iter::once(OsString::from("vep recode")).chain(spelled))
            .unwrap()
    }

    /// `variant_recoder`'s spellings reach the parser: `--id` for the input
    /// string (43), `-i` and `-s`, any case and a unique prefix.
    #[test]
    fn getopt_long_spellings_reach_the_options() {
        let args = parsed(&[
            "--id",
            "rs699",
            "--VCF_STRING",
            "--mane",
            "-s",
            "homo_sapiens",
        ]);
        assert_eq!(args.input_data.as_deref(), Some("rs699"));
        assert!(args.vcf_string && args.mane_select);
        assert_eq!(args.species.as_deref(), Some("homo_sapiens"));
        let args = parsed(&[
            "-i",
            "in.txt",
            "-pretty",
            "--fields=hgvsg,spdi",
            "--pass",
            "x",
        ]);
        assert_eq!(args.input_file.as_deref(), Some("in.txt"));
        assert!(args.pretty);
        assert_eq!(args.fields.as_deref(), Some("hgvsg,spdi"));
        assert_eq!(args.password.as_deref(), Some("x"));
    }

    /// A repeated option is read as Getopt::Long reads `variant_recoder`'s
    /// (36-95): a scalar or switch keeps its last occurrence
    /// (`Getopt/Long.pm:677` assigns each one), and `--transcript_filter`
    /// keeps every occurrence (66).
    #[test]
    fn repeated_options_keep_the_last_value_and_transcript_filters_accumulate() {
        let args = parsed(&[
            "--id",
            "rs699",
            "--id",
            "rs700",
            "--pretty",
            "--pretty",
            "--buffer_size",
            "5",
            "--buffer_size",
            "1",
            "--transcript_filter",
            "stable_id match N[MR]_",
            "--transcript_filter",
            "biotype is protein_coding",
        ]);
        assert_eq!(args.input_data.as_deref(), Some("rs700"));
        assert!(args.pretty);
        assert_eq!(args.buffer_size, Some(1));
        assert_eq!(
            args.transcript_filter,
            ["stable_id match N[MR]_", "biotype is protein_coding"]
        );
        let args = parsed(&[
            "--fields",
            "hgvsg",
            "--fields=hgvsc,spdi",
            "-s",
            "mus",
            "-s",
            "homo_sapiens",
        ]);
        assert_eq!(args.fields.as_deref(), Some("hgvsc,spdi"));
        assert_eq!(args.species.as_deref(), Some("homo_sapiens"));
    }

    /// Every database flag, the LRG and ambiguous-HGVS features and the
    /// alignment file are refused naming the flag; the first given one is named.
    #[test]
    fn database_flags_are_refused_by_name() {
        for argv in [
            vec!["--host", "ensembldb.ensembl.org"],
            vec!["--port", "3306"],
            vec!["-u", "anonymous"],
            vec!["--password", "x"],
            vec!["--registry", "reg.pm"],
            vec!["--db_version", "116"],
            vec!["--genomes"],
            vec!["--is_multispecies", "1"],
            vec!["--lrg"],
            vec!["--ambiguous_hgvs"],
            vec!["--bam", "a.bam"],
        ] {
            let mut full = vec!["--id", "rs699"];
            full.extend(argv.iter().copied());
            let refused = parsed(&full)
                .refused_flag()
                .unwrap_or_else(|| panic!("{} was not refused", full.join(" ")));
            let flag = argv[0].trim_start_matches('-');
            let flag = if flag == "u" { "user" } else { flag };
            assert!(
                refused
                    .message
                    .starts_with(&format!("--{flag} is not supported: it ")),
                "{}",
                refused.message
            );
            assert_eq!(refused.component, "recode");
        }
        assert!(parsed(&["--id", "rs699"]).refused_flag().is_none());
    }
}
