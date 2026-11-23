// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line before clap reads it: Ensembl VEP's configuration sources
//! and its option spellings.
//!
//! Ensembl VEP takes a flag from five sources, the first that names it winning
//! (`Config.pm:693-698` at release/116.3): the command line, the `--config`
//! file, `<dir>/vep.ini`, a `VEP_<FLAG>` environment variable, the defaults.
//! [`prepass`] reads the three sources behind the command line and lays their
//! flags out around it, a flag the command line names dropped from every lower
//! source and a lower source's flag dropped from the ones below it, so clap
//! parses one value per flag and Ensembl VEP's priority holds. Ensembl VEP
//! merges the sources with `||=` (`Config.pm:1101, 1158`), so there a file's
//! value replaces a command-line `0`; here a higher source's value stands
//! whatever it is.
//!
//! Getopt::Long reads an option name in any case, abbreviated to any prefix
//! that names one option, and after one dash as well as two (`Config.pm:62`
//! leaves `ignore_case` and `auto_abbrev` on and `bundling` off); the pre-pass
//! spells each such option as the parser declares it, so `-cache`, `--CACHE`
//! and `--offl` reach `--cache` and `--offline`. A name the parser's own flags leave
//! ambiguous, or one it lacks, is read as Getopt::Long reads it against
//! `@VEP_PARAMS`, so `--js` reaches `--json`, and a one-dash long name such as
//! `-cell_type x` is spelled `--cell_type=x` rather than read as `-c ell_type`.
//! A token no option names is left as written for clap to read or refuse.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::{OsStr, OsString};
use std::fs;
use std::path::{Path, PathBuf};

use clap::{ArgAction, Command};

use crate::args::Args;
use crate::config::{error_table_entry, CliError, ERROR_TABLE};
use crate::option_rules::ALLOW_MULTIPLE;

/// `@VEP_PARAMS` (`Config.pm:73-276`): every option Ensembl VEP 116.3
/// declares, in Getopt::Long's notation (`name|alias=type`) and the source's
/// order. The pre-pass measures its own prefix resolution against this list.
#[doc(hidden)]
pub const VEP_PARAMS: &[&str] = &[
    "help",
    "config=s",
    "input_file|i=s",
    "input_data|id=s",
    "format=s",
    "output_format=s",
    "delimiter=s",
    "no_check_variants_order",
    "species|s=s",
    "registry=s",
    "host=s",
    "port=s",
    "user|u=s",
    "password|pass=s",
    "db_version=i",
    "assembly|a=s",
    "grch37",
    "genomes",
    "refseq",
    "merged",
    "all_refseq",
    "gencode_basic",
    "gencode_primary",
    "is_multispecies=i",
    "transcript_filter=s@",
    "exclude_predicted",
    "minimal",
    "most_severe",
    "summary",
    "pick",
    "pick_allele",
    "per_gene",
    "pick_allele_gene",
    "flag_pick",
    "flag_pick_allele",
    "flag_pick_allele_gene",
    "flag_gencode_primary",
    "pick_order=s",
    "buffer_size=i",
    "failed=i",
    "gp",
    "chr=s",
    "check_ref",
    "lookup_ref",
    "check_existing",
    "check_svs",
    "no_check_alleles",
    "exclude_null_alleles",
    "check_frequency",
    "af",
    "af_1kg",
    "af_gnomade",
    "af_gnomad",
    "af_gnomadg",
    "old_maf",
    "max_af",
    "pubmed",
    "freq_filter=s",
    "freq_freq=f",
    "freq_gt_lt=s",
    "freq_pop=s",
    "filter_common",
    "allow_non_variant",
    "process_ref_homs",
    "individual=s",
    "individual_zyg=s",
    "phased",
    "fork=i",
    "dont_skip",
    "nearest=s",
    "distance=s",
    "clin_sig_allele=i",
    "clinvar_somatic_classification",
    "overlaps",
    "max_sv_size=i",
    "remove_hgvsp_version",
    "verbose|v",
    "quiet|q",
    "no_progress",
    "everything|e",
    "output_file|o=s@",
    "compress_output=s",
    "no_headers",
    "mask_header_cache_path",
    "stats_file|sf=s",
    "stats_text",
    "stats_html",
    "no_stats",
    "warning_file=s",
    "skipped_variants_file=s",
    "force_overwrite|force",
    "terms|t=s",
    "coding_only",
    "canonical",
    "mane",
    "mane_select",
    "tsl",
    "appris",
    "ccds",
    "xref_refseq",
    "uniprot",
    "protein",
    "protein_version",
    "biotype",
    "hgnc",
    "symbol",
    "transcript_version",
    "gene_version",
    "gene_phenotype",
    "mirna",
    "spdi",
    "ga4gh_vrs",
    "hgvs",
    "hgvsg",
    "hgvsg_use_accession",
    "hgvsp_use_prediction",
    "shift_hgvs=i",
    "ambiguous_hgvs",
    "sift=s",
    "polyphen=s",
    "humdiv",
    "condel=s",
    "variant_class",
    "regulatory",
    "cell_type=s",
    "regulatory_gff=s",
    "extended_promoters",
    "convert=s",
    "no_intergenic",
    "vcf",
    "solr",
    "json",
    "tab",
    "vcf_info_field=s",
    "keep_csq",
    "keep_ann",
    "lrg",
    "fields=s",
    "domains",
    "numbers",
    "total_length",
    "allele_number",
    "show_ref_allele",
    "no_escape",
    "ambiguity",
    "var_synonyms",
    "shift_3prime=i",
    "shift_genomic=i",
    "shift_length",
    "uploaded_allele",
    "database",
    "cache:s",
    "cache_version=i",
    "show_cache_info",
    "dir=s",
    "dir_cache=s",
    "full_cache_dir=s",
    "dir_plugins=s",
    "offline",
    "fasta|fa=s",
    "fasta_dir=s",
    "no_fasta",
    "sereal",
    "synonyms=s",
    "bam=s",
    "use_transcript_ref",
    "use_given_ref",
    "custom=s@",
    "tmpdir=s",
    "gff=s",
    "gtf=s",
    "bigwig=s",
    "phyloP=s@",
    "phastCons=s@",
    "ucsc_assembly=s",
    "ucsc_data_root=s",
    "custom_multi_allelic",
    "custom_suppress_filter",
    "plugin=s@",
    "safe",
    "debug",
];

/// The name of the configuration file read from `--dir` (`Config.pm:704`).
#[doc(hidden)]
pub const INI_FILE: &str = "vep.ini";

/// The prefix of an environment variable that sets a flag (`Config.pm:1133`).
#[doc(hidden)]
pub const ENV_PREFIX: &str = "VEP_";

/// The directory `vep.ini` is read from when neither the command line nor the
/// `--config` file names `--dir` (`Config.pm:294, 704`).
#[doc(hidden)]
pub const DEFAULT_DIR: &str = "~/.vep";

/// Flags no source behind the command line can set: `config` names the file
/// being read, and the `vep` script tests `help` before any file is opened
/// (vep:29); `version` is the parser's own.
#[doc(hidden)]
pub const COMMAND_LINE_ONLY: &[&str] = &["config", "help", "version"];

/// The refusal of `-V` on a command line that holds any other token, printed
/// under the `cli` component. `-V` is the version's short, so such a line would
/// print the version and annotate nothing, where Getopt::Long's `ignore_case`
/// (`Config.pm:62`) reads the token as `-v` (`verbose|v`, `:159`) and Ensembl
/// VEP runs.
#[doc(hidden)]
pub const VERSION_SHORT_BESIDE_ARGUMENTS: &str = "-V beside other arguments is not supported: it prints the version here, where Ensembl VEP reads it as -v (--verbose) and runs; write -v for a verbose run, or --version for the version";

/// One option's names: the one it is stored under and every spelling that
/// reaches it (its aliases, and a single-letter short as Getopt::Long reads
/// one, as a name).
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OptionNames {
    pub canonical: String,
    pub names: Vec<String>,
}

/// Whether a written name is one of an option's names: in any case for a long
/// name, exactly for a single letter, which clap reads as a short and tells
/// apart by case.
fn same_name(written: &str, name: &str) -> bool {
    if name.len() == 1 || written.len() == 1 {
        written == name
    } else {
        written.eq_ignore_ascii_case(name)
    }
}

/// The options with a name the written name is a prefix of, in any case.
fn options_prefixed_by<'a>(
    name: &str,
    options: &'a [OptionNames],
) -> impl Iterator<Item = &'a OptionNames> {
    let lower = name.to_ascii_lowercase();
    options.iter().filter(move |o| {
        o.names
            .iter()
            .any(|n| n.to_ascii_lowercase().starts_with(&lower))
    })
}

/// The option a name reaches when read as Getopt::Long reads it: a name that
/// is one of an option's names wins; otherwise the name is a prefix, in any
/// case, of the names of exactly one option; otherwise none.
#[doc(hidden)]
pub fn resolve<'a>(name: &str, options: &'a [OptionNames]) -> Option<&'a OptionNames> {
    if name.is_empty() {
        return None;
    }
    if let Some(exact) = options
        .iter()
        .find(|o| o.names.iter().any(|n| same_name(name, n)))
    {
        return Some(exact);
    }
    let mut prefixed = options_prefixed_by(name, options);
    let first = prefixed.next()?;
    prefixed.next().is_none().then_some(first)
}

/// The value part of the `@VEP_PARAMS` spec naming an option (`=s`, `=i`, `:s`,
/// or empty for a switch); `None` for a name `@VEP_PARAMS` lacks.
fn perl_value_spec(name: &str) -> Option<&'static str> {
    VEP_PARAMS.iter().find_map(|spec| {
        let (names, value_type) = match spec.find(['=', ':']) {
            Some(at) => spec.split_at(at),
            None => (*spec, ""),
        };
        names.split('|').any(|n| n == name).then_some(value_type)
    })
}

/// How Ensembl VEP takes an option's value, from its `@VEP_PARAMS` spec: `=`
/// a required value, `:` an optional one, neither a switch; `None` for a name
/// `@VEP_PARAMS` lacks.
fn perl_arity(name: &str) -> Option<Arity> {
    let value_type = perl_value_spec(name)?;
    let arity = if value_type.is_empty() {
        Arity::Switch
    } else if value_type.starts_with('=') {
        Arity::Required
    } else {
        Arity::Optional
    };
    Some(arity)
}

/// Whether Getopt::Long reads a token asking for the help or the version after
/// the option as that, not as the option's value: the option is declared as a
/// number (`=i`, `=f`), whose value a token that is not a number is pushed back
/// from, or is no option of Ensembl VEP's, which `pass_through` leaves in
/// `@ARGV` with the token after it unread.
fn perl_leaves_a_usage_token_unread(name: &str) -> bool {
    perl_value_spec(name).is_none_or(|spec| spec.ends_with(['i', 'f']))
}

/// Ensembl VEP's options as [`VEP_PARAMS`] declares them.
#[doc(hidden)]
pub fn perl_options() -> Vec<OptionNames> {
    VEP_PARAMS
        .iter()
        .map(|spec| {
            let names = spec
                .split_once(['=', ':'])
                .map_or(*spec, |(names, _value_type)| names);
            let names: Vec<String> = names.split('|').map(str::to_string).collect();
            OptionNames {
                canonical: names[0].clone(),
                names,
            }
        })
        .collect()
}

/// How the parser takes an option's value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Arity {
    /// A switch: `--offline`.
    Switch,
    /// A value may follow: `--cache`, `--cache <dir>`.
    Optional,
    /// A value must follow: `--format <format>`.
    Required,
}

/// An option the parser declares.
#[derive(Debug, Clone)]
struct Declared {
    names: OptionNames,
    arity: Arity,
    /// A repeatable flag that accumulates its values (Ensembl VEP's
    /// `@ALLOW_MULTIPLE`); every source appends to it.
    multiple: bool,
    /// The keys a configuration file or the environment may set the option
    /// under: its long name, and any alias that is an option of Ensembl VEP's
    /// own (`full_cache_dir`, an alias of `json_cache` here, an option of
    /// `@VEP_PARAMS` that Ensembl VEP stores and reads under that name).
    source_keys: Vec<String>,
}

/// The names of Ensembl VEP's options, as `@VEP_PARAMS` stores them.
fn perl_canonicals() -> BTreeSet<String> {
    perl_options().into_iter().map(|o| o.canonical).collect()
}

/// The parser's options, read from the built command.
fn declared_options(command: &Command) -> Vec<Declared> {
    let perl_canonicals = perl_canonicals();
    command
        .get_arguments()
        .filter_map(|arg| {
            let long = arg.get_long()?;
            let mut names = vec![long.to_string()];
            names.extend(
                arg.get_all_aliases()
                    .into_iter()
                    .flatten()
                    .map(str::to_string),
            );
            names.extend(arg.get_short().map(|c| c.to_string()));
            names.extend(
                arg.get_all_short_aliases()
                    .into_iter()
                    .flatten()
                    .map(|c| c.to_string()),
            );
            let arity = match arg.get_num_args() {
                Some(range) if range.takes_values() && range.min_values() == 0 => Arity::Optional,
                Some(range) if range.takes_values() => Arity::Required,
                _ => Arity::Switch,
            };
            let source_keys = names
                .iter()
                .filter(|n| *n == long || perl_canonicals.contains(*n))
                .cloned()
                .collect();
            Some(Declared {
                names: OptionNames {
                    canonical: long.to_string(),
                    names,
                },
                arity,
                multiple: matches!(arg.get_action(), ArgAction::Append)
                    && ALLOW_MULTIPLE.contains(&long),
                source_keys,
            })
        })
        .collect()
}

/// The declared option stored under a long name.
fn declared_option<'a>(declared: &'a [Declared], long: &str) -> Option<&'a Declared> {
    declared.iter().find(|d| d.names.canonical == long)
}

/// The declared option a configuration key names: by its long name, or by an
/// alias that is an Ensembl VEP option of its own; an alias Ensembl VEP also
/// treats as one (`fa`, `force`) names nothing, as its readers store such a key
/// where nothing reads it.
fn declared_for_key<'a>(declared: &'a [Declared], key: &str) -> Option<&'a Declared> {
    declared
        .iter()
        .find(|d| d.source_keys.iter().any(|k| k == key))
}

/// The declared option carrying a resolved name, which may be Ensembl VEP's
/// name for an option the parser declares under another (`full_cache_dir` for
/// `json_cache`).
fn declared_by_name<'a>(declared: &'a [Declared], name: &str) -> Option<&'a Declared> {
    declared
        .iter()
        .find(|d| d.names.names.iter().any(|n| same_name(name, n)))
}

/// The parser's options, each with the names that reach it.
#[doc(hidden)]
pub fn parser_options(command: &Command) -> Vec<OptionNames> {
    declared_options(command)
        .into_iter()
        .map(|d| d.names)
        .collect()
}

/// A companion tool's command line spelled as its parser declares it, as
/// Getopt::Long reads `filter_vep`'s and `variant_recoder`'s: an option
/// written after one dash or two, in any case, or as a prefix of one option's
/// names, becomes `--<long>` with its attached value; the token after an
/// option that takes a value is that value however it is spelled; `--` ends
/// the options. A name of two or more letters that is a prefix of several
/// options is written with two dashes, so clap refuses it rather than reading
/// its letters as short forms (`-li`). Any other name no option has is left as
/// written for clap to refuse.
#[doc(hidden)]
pub fn spell_as_declared(argv: &[OsString], command: &Command) -> Vec<OsString> {
    let options = parser_options(command);
    let takes_value = |long: &str| {
        command
            .get_arguments()
            .find(|arg| arg.get_long() == Some(long))
            .is_some_and(|arg| arg.get_num_args().is_some_and(|n| n.takes_values()))
    };
    let mut out = Vec::with_capacity(argv.len());
    let mut value_next = false;
    let mut after_separator = false;
    for token in argv {
        if value_next || after_separator {
            value_next = false;
            out.push(token.clone());
            continue;
        }
        let Some(text) = token.to_str() else {
            out.push(token.clone());
            continue;
        };
        if text == "--" {
            after_separator = true;
            out.push(token.clone());
            continue;
        }
        let Some((name, attached)) = option_parts(text) else {
            out.push(token.clone());
            continue;
        };
        match resolve(name, &options) {
            Some(option) => {
                let long = option.canonical.as_str();
                out.push(OsString::from(match attached {
                    Some(value) => format!("--{long}={value}"),
                    None => format!("--{long}"),
                }));
                value_next = attached.is_none() && takes_value(long);
            }
            // Getopt::Long calls a prefix of several options ambiguous whatever
            // its dashes; written with two, clap refuses the name instead of
            // reading `-li` as `-l -i`.
            None if name.chars().nth(1).is_some()
                && options_prefixed_by(name, &options).next().is_some() =>
            {
                let two_dashes = Token::option(name, attached.map(OsString::from));
                out.push(two_dashes.into_os_string());
            }
            None => out.push(token.clone()),
        }
    }
    out
}

/// The names the binary knows: every declared option, and every flag the error
/// table refuses before parsing, which [`crate::config::parse_command_line`]
/// finds by its double-dash spelling.
fn known_options(declared: &[Declared]) -> Vec<OptionNames> {
    let mut known: Vec<OptionNames> = declared.iter().map(|d| d.names.clone()).collect();
    known.extend(
        ERROR_TABLE
            .iter()
            .filter(|entry| !entry.declared)
            .map(|entry| OptionNames {
                canonical: entry.flag.to_string(),
                names: vec![entry.flag.to_string()],
            }),
    );
    known
}

/// The option a written name reaches: one the binary knows, by [`resolve`];
/// failing that, the one Getopt::Long reaches against `@VEP_PARAMS`, so a
/// prefix a flag vep-rs adds makes ambiguous (`js`: `json` and `json_cache`)
/// reads as Ensembl VEP reads it, and an Ensembl VEP option the binary lacks
/// is named in full for the parser to refuse.
#[doc(hidden)]
pub fn resolve_spelling<'a>(
    name: &str,
    known: &'a [OptionNames],
    perl: &'a [OptionNames],
) -> Option<&'a OptionNames> {
    resolve(name, known).or_else(|| resolve(name, perl))
}

/// A spelling Getopt::Long resolves to one Ensembl VEP option and the pre-pass
/// does not: the option's own name or alias, or a prefix of one, that resolves
/// to another option here or to none.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnresolvedSpelling {
    /// The spelling as Getopt::Long takes it, without its dashes.
    pub spelling: String,
    /// The option it reaches in Ensembl VEP.
    pub perl: String,
    /// The parser's options it is a prefix of, sorted; empty when none is.
    pub candidates: Vec<String>,
}

/// Every spelling of an Ensembl VEP option the parser declares that
/// Getopt::Long resolves and [`resolve_spelling`] does not, sorted by spelling;
/// empty while every such spelling reaches the same option.
#[doc(hidden)]
pub fn unresolved_perl_spellings(command: &Command) -> Vec<UnresolvedSpelling> {
    let perl = perl_options();
    let declared = declared_options(command);
    let known = known_options(&declared);
    let parser = parser_options(command);
    let mut out = Vec::new();
    for option in &perl {
        if !parser.iter().any(|p| p.canonical == option.canonical) {
            continue;
        }
        for name in &option.names {
            for end in 1..=name.len() {
                let spelling = &name[..end];
                let getopt_long_reads_it =
                    resolve(spelling, &perl).is_some_and(|o| o.canonical == option.canonical);
                let pre_pass_reads_it = resolve_spelling(spelling, &known, &perl)
                    .is_some_and(|o| o.canonical == option.canonical);
                if !getopt_long_reads_it || pre_pass_reads_it {
                    continue;
                }
                let mut candidates: Vec<String> = options_prefixed_by(spelling, &parser)
                    .map(|p| p.canonical.clone())
                    .collect();
                candidates.sort();
                candidates.dedup();
                out.push(UnresolvedSpelling {
                    spelling: spelling.to_string(),
                    perl: option.canonical.clone(),
                    candidates,
                });
            }
        }
    }
    out.sort_by(|a, b| a.spelling.cmp(&b.spelling));
    out.dedup();
    out
}

/// Perl's truth of a configuration value: anything but the empty string and `0`.
pub(crate) fn truthy(value: &str) -> bool {
    !(value.is_empty() || value == "0")
}

/// The integer Getopt::Long stores for an `=i` value: an optional sign, then
/// digits with `_` allowed before and among them and dropped, read as `0+$arg`,
/// so `+1`, `01` and `_1` are 1 and `-0` and `0_0` are 0. One line feed closing
/// the text is passed over, as the `$` of Getopt::Long's pattern passes it. A
/// magnitude beyond `i64::MAX` reads as that bound, non-zero as the float Perl
/// stores; `None` for a text Getopt::Long refuses as the option's value.
pub(crate) fn getopt_integer(text: &str) -> Option<i64> {
    let text = text.strip_suffix('\n').unwrap_or(text);
    let (negative, unsigned) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text.strip_prefix('+').unwrap_or(text)),
    };
    let digits = unsigned.trim_start_matches('_');
    let is_integer = digits.starts_with(|c: char| c.is_ascii_digit())
        && digits.bytes().all(|b| b.is_ascii_digit() || b == b'_');
    if !is_integer {
        return None;
    }
    let magnitude = digits
        .bytes()
        .filter(u8::is_ascii_digit)
        .fold(0i64, |number, digit| {
            number
                .saturating_mul(10)
                .saturating_add(i64::from(digit - b'0'))
        });
    Some(if negative { -magnitude } else { magnitude })
}

/// The `=i` options Ensembl VEP reads for truth alone (`shift_3prime=i`,
/// `shift_genomic=i`, `Config.pm:231-232`). The command line's value is the
/// number Getopt::Long stores, while a value from a file or the environment
/// stays the text its source holds (`Config.pm:1101, 1158`), true unless it is
/// `0`, so `00` is off on the command line and on from a file.
const READ_FOR_TRUTH: &[&str] = &["shift_3prime", "shift_genomic"];

/// Whether a key names a flag the error table refuses before parsing, one the
/// parser does not declare ([`known_options`] lists them).
fn refused_before_parsing(key: &str) -> bool {
    error_table_entry(key).is_some_and(|entry| !entry.declared)
}

/// Where a source's key is stored, and whether it appends: a declared option's
/// long name and its `multiple` ([`declared_for_key`]); for an Ensembl VEP
/// option the parser lacks, or a flag the error table refuses before parsing,
/// the flag's own name, so the parser refuses it as it refuses the flag on the
/// command line, appending when `@ALLOW_MULTIPLE` lists it; `None` for a key
/// none of them names.
fn key_storage(
    declared: &[Declared],
    perl: &BTreeSet<String>,
    key: &str,
) -> Option<(String, bool)> {
    match declared_for_key(declared, key) {
        Some(option) => Some((option.names.canonical.clone(), option.multiple)),
        None if perl.contains(key) || refused_before_parsing(key) => {
            Some((key.to_string(), ALLOW_MULTIPLE.contains(&key)))
        }
        None => None,
    }
}

/// The flags one source contributes, by long name: a scalar flag holds one
/// value, a repeatable one every value in order.
#[derive(Debug, Default)]
struct Settings {
    scalars: BTreeMap<String, String>,
    lists: BTreeMap<String, Vec<String>>,
    /// The key a flag was set under where it is Ensembl VEP's own name for an
    /// option declared here as an alias (`full_cache_dir` for `json_cache`), the
    /// name its token is written under; absent for a flag set under its long name.
    written: BTreeMap<String, String>,
}

impl Settings {
    /// Records one `key value` pair as Ensembl VEP's readers do
    /// (`Config.pm:1094-1102`, `:1150-1158`): a key that only the command line
    /// may set, or that [`key_storage`] stores nowhere, is dropped; a repeatable
    /// key appends; any other key keeps the source's first value that is not `0`
    /// or empty (`||=`).
    fn record(&mut self, declared: &[Declared], perl: &BTreeSet<String>, key: &str, value: &str) {
        if COMMAND_LINE_ONLY.contains(&key) {
            return;
        }
        let Some((long, multiple)) = key_storage(declared, perl, key) else {
            return;
        };
        let written = (key != long && perl.contains(key)).then(|| key.to_string());
        if multiple {
            self.lists
                .entry(long.clone())
                .or_default()
                .push(value.to_string());
        } else if !self
            .scalars
            .get(&long)
            .is_some_and(|current| truthy(current))
        {
            self.scalars.insert(long.clone(), value.to_string());
        } else {
            return;
        }
        match written {
            Some(written) => self.written.insert(long, written),
            None => self.written.remove(&long),
        };
    }

    /// The name a stored flag's token is written under ([`Settings::written`]).
    fn written_name<'a>(&'a self, long: &'a str) -> &'a str {
        self.written.get(long).map_or(long, String::as_str)
    }
}

/// Whether a character is a blank to Perl's `\s` on the bytes of a configuration
/// line: the ASCII blanks alone.
fn is_blank(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0B' | '\x0C')
}

/// The words of a configuration line as `read_config_from_file` splits them
/// (`Config.pm:1082-1092`): a run of blanks followed by an odd number of double
/// quotes on the line is one blank inside its word; every other run separates
/// words; then every `"` and `'` is removed from the words after the first. A
/// line starting with a blank has an empty first word, which no option is
/// named by.
fn split_config_line(line: &str) -> Vec<String> {
    let total_quotes = line.matches('"').count();
    let mut words = vec![String::new()];
    let mut quotes_seen = 0usize;
    let mut in_blank = false;
    for c in line.chars() {
        if is_blank(c) {
            if !in_blank {
                if (total_quotes - quotes_seen) % 2 == 1 {
                    words.last_mut().expect("one word").push(' ');
                } else {
                    words.push(String::new());
                }
            }
            in_blank = true;
            continue;
        }
        in_blank = false;
        if c == '"' {
            quotes_seen += 1;
        }
        words.last_mut().expect("one word").push(c);
    }
    if words.len() > 1 && words.last().is_some_and(String::is_empty) {
        words.pop();
    }
    for word in words.iter_mut().skip(1) {
        word.retain(|c| c != '"' && c != '\'');
    }
    words
}

/// Reads one configuration file as `read_config_from_file` does
/// (`Config.pm:1068-1110`): a path Perl reads as false (`0`, the empty value)
/// is refused before any open, shown as the empty string (`:1073`); a
/// directory opens and yields no line, so nothing is read (`:1075-1077`); a
/// line starting with `#`, or blank, is skipped; the first word, less one
/// leading dash, is the key; the value of a repeatable key is every following
/// word joined by one blank, of any other key the second word alone; a key
/// with no value is dropped, as a key that names no option of Ensembl VEP's or
/// the parser's is. The file is read as bytes, as Ensembl VEP reads it; a byte that
/// is not UTF-8 becomes the replacement character in a value and is skipped in
/// a comment.
fn read_config_file(
    path: &Path,
    shown: &str,
    declared: &[Declared],
    perl: &BTreeSet<String>,
) -> Result<Settings, CliError> {
    if !truthy(shown) {
        return Err(CliError::new("config", "Could not open config file \"\""));
    }
    if path.is_dir() {
        return Ok(Settings::default());
    }
    let bytes = fs::read(path)
        .map_err(|_| CliError::new("config", format!("Could not open config file \"{shown}\"")))?;
    let text = String::from_utf8_lossy(&bytes);
    let mut settings = Settings::default();
    for line in text.lines() {
        if line.starts_with('#') || line.chars().all(is_blank) {
            continue;
        }
        let words = split_config_line(line);
        let key = words[0].strip_prefix('-').unwrap_or(&words[0]);
        let values = &words[1..];
        if values.is_empty() {
            continue;
        }
        let multiple = key_storage(declared, perl, key).is_some_and(|(_, multiple)| multiple);
        let value = if multiple {
            values.join(" ")
        } else {
            values[0].clone()
        };
        settings.record(declared, perl, key, &value);
    }
    Ok(settings)
}

/// The flags of the environment, as `read_config_from_environment` reads them
/// (`Config.pm:1127-1168`): every `VEP_<NAME>` variable with a non-empty value
/// sets the flag `<name>`, lower-cased; the variables are taken in name order,
/// where Ensembl VEP takes them in hash order. The lower-cased name is matched
/// exactly, so no variable names `phyloP` or `phastCons`, the two keys spelled
/// with a capital letter: Ensembl VEP stores such a variable under the
/// lower-cased key (`:1140, 1158`), which its option sets never read
/// (`:553, 560`).
fn read_environment(
    env: &BTreeMap<String, String>,
    declared: &[Declared],
    perl: &BTreeSet<String>,
) -> Settings {
    let mut settings = Settings::default();
    for (var, value) in env {
        let Some(name) = var.strip_prefix(ENV_PREFIX) else {
            continue;
        };
        if value.is_empty() {
            continue;
        }
        settings.record(declared, perl, &name.to_ascii_lowercase(), value);
    }
    settings
}

/// A token of the command line after spelling.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    /// An option the parser declares, by its long name, with its value: the one
    /// attached by `=`, or the token after an option that requires one.
    Option {
        long: String,
        /// The name the token is written under for clap: the long name, or
        /// Ensembl VEP's own name for an option declared here as an alias of
        /// another (`full_cache_dir` for `json_cache`), which is the key Ensembl
        /// VEP stores the flag as and prints in its command-line header
        /// (`Config.pm:734-745`).
        written: String,
        value: Option<OsString>,
        /// Whether the token is the version asked by its one-letter short, `-V`,
        /// the token Getopt::Long's `ignore_case` reads as `-v`
        /// (`Config.pm:62, 159`).
        version_short: bool,
    },
    /// A token left as written: a value, a short cluster, the `--` separator, an
    /// option the parser does not declare, or anything after `--`.
    Verbatim(OsString),
}

impl Token {
    /// An option written under its long name.
    fn option(long: &str, value: Option<OsString>) -> Token {
        Token::Option {
            long: long.to_string(),
            written: long.to_string(),
            value,
            version_short: false,
        }
    }

    /// The token as clap reads it; a value is attached with `=`, so a value
    /// starting with a dash stays a value, as Getopt::Long takes it.
    fn into_os_string(self) -> OsString {
        match self {
            Token::Option { written, value, .. } => {
                let mut out = OsString::from(format!("--{written}"));
                if let Some(value) = value {
                    out.push("=");
                    out.push(value);
                }
                out
            }
            Token::Verbatim(s) => s,
        }
    }
}

/// Ensembl VEP's own name for an option the parser declares under another: the
/// `@VEP_PARAMS` option the written name resolves to, when it is one of the
/// parser's names for the option but not its long name (`full_cache_dir`, an
/// alias of `json_cache`); `None` when the long name is the key both store.
fn perl_alias_key(name: &str, carrier: &Declared, perl: &[OptionNames]) -> Option<String> {
    let key = &resolve(name, perl)?.canonical;
    let is_alias = key != &carrier.names.canonical && carrier.names.names.contains(key);
    is_alias.then(|| key.clone())
}

/// Whether a token is a bare option asking for the help or the version, in any
/// spelling the pre-pass reads (`--help`, `-h`, `-help`, `--HELP`, `--hel`,
/// `--version`, `-V`).
fn asks_usage(token: &OsStr, known: &[OptionNames], perl: &[OptionNames]) -> bool {
    let Some((name, None)) = token.to_str().and_then(option_parts) else {
        return false;
    };
    resolve_spelling(name, known, perl)
        .is_some_and(|option| option.canonical == "help" || option.canonical == "version")
}

/// The option name and attached value a token carries when it is written as
/// one: `--name`, `-name`, `--name=value`, `-name=value`; `--` alone, a lone
/// dash and anything else is not an option.
fn option_parts(token: &str) -> Option<(&str, Option<&str>)> {
    let body = token
        .strip_prefix("--")
        .or_else(|| token.strip_prefix('-'))?;
    if body.is_empty() {
        return None;
    }
    Some(match body.split_once('=') {
        Some((name, value)) => (name, Some(value)),
        None => (body, None),
    })
}

/// A short flag with its value attached, as clap reads one (`-ofile.txt`): the
/// first letter of the token's body is a declared short that requires a value,
/// and the rest of the body is the value. Getopt::Long refuses the spelling, so
/// a token that names an Ensembl VEP option is never read this way.
fn glued_short<'a, 'n>(body: &'n str, declared: &'a [Declared]) -> Option<(&'a Declared, &'n str)> {
    let mut chars = body.chars();
    let short = chars.next()?.to_string();
    let rest = chars.as_str();
    if rest.is_empty() {
        return None;
    }
    declared
        .iter()
        .find(|d| d.arity == Arity::Required && d.names.names.contains(&short))
        .map(|d| (d, rest))
}

/// Spells the command line's own tokens: an option that [`resolve_spelling`]
/// reaches is written as `--<long>`, with its attached value, under Ensembl
/// VEP's own name where that is an alias here ([`perl_alias_key`]); the token
/// after an option that requires a value is that value, however it is spelled
/// (`--input_file -cache` names the file `-cache`, as Getopt::Long reads it),
/// except that an option Ensembl VEP declares as a number, or does not declare,
/// leaves a token asking for the help or the version unread, as Getopt::Long
/// pushes back a value that is not a number and passes an unknown option
/// through (`--fork --help` and `--pipeline_depth --help` ask for the help); an
/// option Ensembl VEP declares and the parser does not takes its value as
/// `@VEP_PARAMS` declares it (`--bigwig --help` names the file `--help`), so the
/// parser refuses it by name; a name that is a prefix of several options is
/// written with two dashes, so clap refuses it; a short flag with its value
/// attached is written as the long form; `--` ends the options. Anything else,
/// a cluster of switches or a name no option has, is left for clap. The
/// version written as its one-letter short is marked as that.
fn spell(options: &[OsString], declared: &[Declared]) -> Vec<Token> {
    let known = known_options(declared);
    let perl = perl_options();
    let mut tokens = Vec::with_capacity(options.len());
    let mut after_separator = false;
    let mut iter = options.iter().peekable();
    while let Some(token) = iter.next() {
        let Some(text) = token.to_str().filter(|_| !after_separator) else {
            tokens.push(Token::Verbatim(token.clone()));
            continue;
        };
        if text == "--" {
            after_separator = true;
            tokens.push(Token::Verbatim(token.clone()));
            continue;
        }
        let Some((name, attached)) = option_parts(text) else {
            tokens.push(Token::Verbatim(token.clone()));
            continue;
        };
        if let Some(option) = resolve_spelling(name, &known, &perl) {
            let carrier = declared_by_name(declared, &option.canonical);
            let (long, arity) = match carrier {
                Some(carrier) => (carrier.names.canonical.clone(), carrier.arity),
                None => (
                    option.canonical.clone(),
                    perl_arity(&option.canonical).unwrap_or(Arity::Switch),
                ),
            };
            let written = carrier
                .and_then(|carrier| perl_alias_key(name, carrier, &perl))
                .unwrap_or_else(|| long.clone());
            let usage_token_pushed_back = perl_leaves_a_usage_token_unread(&written)
                && iter
                    .peek()
                    .is_some_and(|next| asks_usage(next, &known, &perl));
            let value = match attached {
                Some(attached) => Some(OsString::from(attached)),
                None if arity == Arity::Required && !usage_token_pushed_back => {
                    iter.next().cloned()
                }
                None => None,
            };
            let version_short = long == "version" && name.len() == 1;
            tokens.push(Token::Option {
                long,
                written,
                value,
                version_short,
            });
            continue;
        }
        // A name that is a prefix of several options is ambiguous to Getopt::Long
        // whatever its dashes; written with two, clap refuses it instead of
        // reading `-sh` as `-s h`.
        let prefixes_an_option = options_prefixed_by(name, &known).next().is_some()
            || options_prefixed_by(name, &perl).next().is_some();
        if prefixes_an_option {
            let two_dashes = Token::option(name, attached.map(OsString::from)).into_os_string();
            tokens.push(Token::Verbatim(two_dashes));
            continue;
        }
        if let Some(body) = text.strip_prefix('-').filter(|b| !b.starts_with('-')) {
            if let Some((option, value)) = glued_short(body, declared) {
                tokens.push(Token::option(
                    &option.names.canonical,
                    Some(OsString::from(value)),
                ));
                continue;
            }
        }
        tokens.push(Token::Verbatim(token.clone()));
    }
    tokens
}

/// Whether the command line asks for the version by `-V` beside another token:
/// the line holds more than one token, and every token of it that asks for the
/// help or the version is the version's one-letter short, so none spells either
/// out, which would win over every other token. A `-V` that is a flag's value,
/// stands behind `--` or is a letter of a cluster is no such token.
fn version_short_beside_arguments(tokens: &[Token]) -> bool {
    let mut asked_by_the_short = tokens.iter().filter_map(|t| match t {
        Token::Option {
            long,
            value: None,
            version_short,
            ..
        } if long == "help" || long == "version" => Some(*version_short),
        _ => None,
    });
    tokens.len() > 1
        && asked_by_the_short.next() == Some(true)
        && asked_by_the_short.all(|short| short)
}

/// The value the command line gives an option: the last occurrence's.
fn command_line_value(tokens: &[Token], long: &str) -> Option<OsString> {
    tokens
        .iter()
        .rev()
        .find_map(|t| match t {
            Token::Option {
                long: name, value, ..
            } if name == long => Some(value),
            _ => None,
        })
        .and_then(|value| value.clone())
}

/// How a source's key takes its value: the declared option's arity, else the
/// `@VEP_PARAMS` spec's for an Ensembl VEP option the parser lacks, else a
/// switch's.
fn key_arity(declared: &[Declared], key: &str) -> Arity {
    declared_option(declared, key)
        .map(|d| d.arity)
        .or_else(|| perl_arity(key))
        .unwrap_or(Arity::Switch)
}

/// The tokens a source's settings add: a switch whose value is true as the flag
/// alone, a flag taking a value as `--<long>=<value>`, a repeatable flag once
/// per value, each under the name it was set by ([`Settings::written_name`]).
/// A switch whose value is `0` or empty, or `cache` set to `0`, adds nothing:
/// Ensembl VEP stores the value and reads it as false. `cache` set to the empty
/// string adds the flag alone: `Config.pm:716` compares the stored value with `0`
/// as strings, which only `0` fails, and `:718` sets it to `1`. A scalar of an
/// Ensembl VEP option the parser lacks is added as its `@VEP_PARAMS` spec takes
/// it, a switch's or a value's, for the parser to refuse; one of a flag the
/// error table refuses before parsing is added as a switch, for the scan ahead
/// of the parser to refuse. A [`READ_FOR_TRUTH`] value that is true as text and
/// zero as a number is added as `1`, the truth Ensembl VEP reads from it.
fn tokens_of(settings: &Settings, declared: &[Declared]) -> (Vec<OsString>, Vec<OsString>) {
    let mut scalars = Vec::new();
    for (key, value) in &settings.scalars {
        let name = settings.written_name(key);
        let token = match (key_arity(declared, key), value.as_str()) {
            (Arity::Switch, v) if !truthy(v) => continue,
            (Arity::Optional, "0") => continue,
            (Arity::Switch, _) | (Arity::Optional, "") => format!("--{name}"),
            (Arity::Required, v)
                if READ_FOR_TRUTH.contains(&key.as_str())
                    && truthy(v)
                    && getopt_integer(v) == Some(0) =>
            {
                format!("--{name}=1")
            }
            (Arity::Optional | Arity::Required, v) => format!("--{name}={v}"),
        };
        scalars.push(token.into());
    }
    let lists = settings
        .lists
        .iter()
        .flat_map(|(key, values)| {
            let name = settings.written_name(key);
            values.iter().map(move |v| format!("--{name}={v}").into())
        })
        .collect();
    (scalars, lists)
}

/// The switches Ensembl VEP reads by `defined()` rather than by truth, so a
/// source's `0` or empty value sets them off but present: `var_synonyms`, which
/// `Runner.pm:117` tests with `defined` before switching `check_existing` on.
pub(crate) const DEFINED_SWITCHES: &[&str] = &["var_synonyms"];

/// The [`DEFINED_SWITCHES`] a source set to a false value, each with that value
/// (`0` or the empty string a quoted value leaves), which [`tokens_of`] adds
/// nothing for; `Config.pm:1101` and `:1158` store the value under the key.
fn defined_off_switches(sources: &[Settings], declared: &[Declared]) -> Vec<(String, String)> {
    DEFINED_SWITCHES
        .iter()
        .filter(|flag| key_arity(declared, flag) == Arity::Switch)
        .filter_map(|flag| {
            sources.iter().find_map(|source| {
                source
                    .scalars
                    .get(*flag)
                    .filter(|value| !truthy(value))
                    .map(|value| (flag.to_string(), value.clone()))
            })
        })
        .collect()
}

/// A path with a leading `~` expanded to `home`, the home directory, which
/// Ensembl VEP leaves to the shell; with no `HOME` (`home` empty), `~/.vep` is
/// `/.vep`, as `$ENV{HOME}.'/.vep'` is (`Config.pm:294`), and `~` alone is `/`.
pub(crate) fn expand_tilde_text(path: &str, home: &str) -> String {
    if path == "~" {
        return if home.is_empty() {
            "/".to_string()
        } else {
            home.to_string()
        };
    }
    match path.strip_prefix("~/") {
        Some(rest) => format!("{home}/{rest}"),
        None => path.to_string(),
    }
}

/// [`expand_tilde_text`] over an argument; a path that is not UTF-8 is taken as
/// written.
fn expand_tilde(path: &OsStr, home: &str) -> PathBuf {
    match path.to_str() {
        Some(text) => PathBuf::from(expand_tilde_text(text, home)),
        None => PathBuf::from(path),
    }
}

/// The environment as `%ENV` holds it, read as bytes (Config.pm:1137 skips the
/// empty value alone); a byte that is not UTF-8 becomes the replacement
/// character, as in a configuration file.
fn environment_map<I>(env: I) -> BTreeMap<String, String>
where
    I: IntoIterator<Item = (OsString, OsString)>,
{
    env.into_iter()
        .map(|(k, v)| {
            (
                k.to_string_lossy().into_owned(),
                v.to_string_lossy().into_owned(),
            )
        })
        .collect()
}

/// The three sources behind a command line, read as `Config->new` reads them
/// (`Config.pm:701-708`): the `--config` file when `config` names one, then
/// `<dir>/vep.ini` when it is a file, `<dir>` being `dir` when it reads as true
/// (`$config->{dir} || $DEFAULTS{dir}`, `:704`), else the `--config` file's
/// `dir`, else `~/.vep`, then the `VEP_*` variables of `env` (`:707`). A flag
/// `named` by the command line is dropped from every source, and a flag a
/// source names from the ones below it, so one value per flag remains and
/// Ensembl VEP's priority holds. `HOME` is read from `env`.
///
/// # Errors
///
/// A `--config` file, or a `vep.ini` that is a file, that cannot be read
/// ([`read_config_file`]).
fn sources_behind(
    config: Option<&OsStr>,
    dir: Option<&OsStr>,
    named: &BTreeSet<String>,
    env: &BTreeMap<String, String>,
    declared: &[Declared],
    perl: &BTreeSet<String>,
) -> Result<[Settings; 3], CliError> {
    let home = env.get("HOME").map_or("", String::as_str);
    let config_file = match config {
        Some(shown) => read_config_file(
            &expand_tilde(shown, home),
            &shown.to_string_lossy(),
            declared,
            perl,
        )?,
        None => Settings::default(),
    };
    let ini_dir = dir
        .filter(|d| truthy(&d.to_string_lossy()))
        .map(OsStr::to_os_string)
        .or_else(|| {
            config_file
                .scalars
                .get("dir")
                .filter(|d| truthy(d))
                .map(OsString::from)
        })
        .unwrap_or_else(|| OsString::from(DEFAULT_DIR));
    let ini_path = expand_tilde(&ini_dir, home).join(INI_FILE);
    // In Ensembl VEP, Perl's `-e` then `open` on a directory reads nothing from it.
    let ini = if ini_path.is_file() {
        read_config_file(&ini_path, &ini_path.display().to_string(), declared, perl)?
    } else {
        Settings::default()
    };
    let environment = read_environment(env, declared, perl);

    let mut named_by_a_higher_source = named.clone();
    let mut sources = [config_file, ini, environment];
    for source in sources.iter_mut() {
        source
            .scalars
            .retain(|key, _| !named_by_a_higher_source.contains(key));
        named_by_a_higher_source.extend(source.scalars.keys().cloned());
    }
    Ok(sources)
}

/// The command line clap parses, laid out as `Config->new` leaves the flags
/// (`Config.pm:701-708`): the sources' scalar flags, lowest source first, then
/// the command line's own tokens, then the sources' repeatable flags in Ensembl
/// VEP's order, command line first (`:757`, `:1097`, `:1152`), then the tokens
/// from a `--` separator on (`rest`).
fn lay_out(
    [config_file, ini, environment]: &[Settings; 3],
    own: &[OsString],
    rest: &[OsString],
    declared: &[Declared],
) -> Vec<OsString> {
    let (config_scalars, config_lists) = tokens_of(config_file, declared);
    let (ini_scalars, ini_lists) = tokens_of(ini, declared);
    let (env_scalars, env_lists) = tokens_of(environment, declared);
    let mut argv = env_scalars;
    argv.extend(ini_scalars);
    argv.extend(config_scalars);
    argv.extend(own.iter().cloned());
    argv.extend(config_lists);
    argv.extend(ini_lists);
    argv.extend(env_lists);
    argv.extend(rest.iter().cloned());
    argv
}

/// A companion tool's command line, spelled by [`spell_as_declared`], with the
/// sources behind it read as `Config->new` reads them for `variant_recoder`
/// (`Config.pm:701-708`, `:1068-1110`, `sources_behind`): the `--config` file
/// it names, `<dir>/vep.ini` under the file's `dir` or `~/.vep`, then the
/// `VEP_*` variables of `env`, laid out around the command line by
/// `lay_out`, the command line keeping its own value of a flag a source
/// names (`:1101`, `:1158`). A key a source names that `command` has no option
/// for is dropped, as `is_valid_param` drops a key `@VEP_PARAMS` lacks
/// (`:1096`, `:1100`), a `config` key inside a file is dropped too, and so is
/// every key of `preset`: the options the tool sets before `Config->new` runs,
/// which Perl's `||=` leaves as set, and the ones it reads before the sources
/// are open or does not know, so that a source's value never reaches them. A
/// leading `~` in the file's path is expanded to the `HOME` of `env`. A command
/// line naming `--help` comes back as spelled, no source read
/// (`variant_recoder:97`).
///
/// # Errors
///
/// A `--config` file, or a `vep.ini` that is a file, that cannot be read:
/// `ERROR: [config] Could not open config file "<path>"` (`Config.pm:1075`).
pub(crate) fn with_config_file<I>(
    command: &Command,
    spelled: &[OsString],
    env: I,
    preset: &[&str],
) -> Result<Vec<OsString>, CliError>
where
    I: IntoIterator<Item = (OsString, OsString)>,
{
    let separator = spelled
        .iter()
        .position(|t| t == "--")
        .unwrap_or(spelled.len());
    let (own, rest) = spelled.split_at(separator);
    let declared = declared_options(command);
    let mut config_path = None;
    let mut named: BTreeSet<String> = BTreeSet::new();
    let mut tokens = own.iter();
    while let Some(token) = tokens.next() {
        // An option `spell_as_declared` resolved is written `--<long>` or
        // `--<long>=<value>`; every other token is a value or left for clap.
        let Some((name, attached)) = token
            .to_str()
            .filter(|text| text.starts_with("--"))
            .and_then(option_parts)
        else {
            continue;
        };
        let takes_value =
            declared_option(&declared, name).is_some_and(|option| option.arity != Arity::Switch);
        let value = match attached {
            Some(value) => Some(OsString::from(value)),
            None if takes_value => tokens.next().cloned(),
            None => None,
        };
        if name == "config" {
            config_path = value;
        }
        named.insert(name.to_string());
    }
    // `variant_recoder:97` prints its usage and exits before `Config->new`
    // reads any file, so a `--help` leaves every source unread.
    if named.contains("help") {
        return Ok(spelled.to_vec());
    }
    // The sources are read with Ensembl VEP's own names, so a file's `dir`
    // places `vep.ini` as `:704` places it; then every key `command` has no
    // option for is dropped, the flags the error table refuses before parsing
    // included, since the companion tool runs no such scan.
    let sources = sources_behind(
        config_path.as_deref(),
        None,
        &named,
        &environment_map(env),
        &declared,
        &perl_canonicals(),
    )?;
    let settable = |key: &String| {
        declared_for_key(&declared, key).is_some() && !preset.contains(&key.as_str())
    };
    let sources = sources.map(|mut source| {
        source.scalars.retain(|key, _| settable(key));
        source.lists.retain(|key, _| settable(key));
        source
    });
    Ok(lay_out(&sources, own, rest, &declared))
}

/// The command line clap parses, from the one given and the sources behind it
/// (`Prepassed::argv` of `prepass_with_sources`).
///
/// The first element is the program name. A command line asking for the help or
/// the version (`--help`, `-h`, `--version`, `-V`, in any spelling) is returned
/// spelled and alone, as the `vep` script answers it before any file is read
/// (vep:29). Otherwise the `--config` file the command line names is read
/// (`Config.pm:701`), then `<dir>/vep.ini` when it is a file, `<dir>` being the
/// command line's `--dir`, else the `--config` file's, else `~/.vep`
/// (`:704-705`), then the `VEP_*` variables of `env` (`:707-708`); a flag a
/// higher source names is dropped from every lower one. The result is the
/// sources' scalar flags, lowest source first, then the command line, then
/// the sources' repeatable flags in Ensembl VEP's order, command line first;
/// the repeatable flags of the sources come before a `--` separator when the
/// command line has one. `HOME` is read from `env`.
///
/// # Errors
///
/// A `--config` file, or a `vep.ini` that is a file, that cannot be read:
/// `ERROR: [config] Could not open config file "<path>"` (`Config.pm:1075`).
/// A `-V` beside any other token, unless the command line spells `--help` or
/// `--version` out: `ERROR: [cli]` with `VERSION_SHORT_BESIDE_ARGUMENTS`,
/// before any source is read.
pub fn prepass<I>(argv: Vec<OsString>, env: I) -> Result<Vec<OsString>, CliError>
where
    I: IntoIterator<Item = (OsString, OsString)>,
{
    prepass_with_sources(argv, env).map(|prepassed| prepassed.argv)
}

/// A command line read through the pre-pass.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Prepassed {
    /// The command line clap parses ([`prepass`]).
    pub argv: Vec<OsString>,
    /// The [`DEFINED_SWITCHES`] a source behind the command line set to `0` or
    /// the empty string, each with that value: off, so no token of `argv` names
    /// them, but present in Ensembl VEP's configuration, where `Config.pm:1101`
    /// and `:1158` store the value.
    pub defined_off: Vec<(String, String)>,
}

/// [`prepass`] with the flags the sources set to a defined false value.
///
/// # Errors
///
/// As [`prepass`].
#[doc(hidden)]
pub fn prepass_with_sources<I>(argv: Vec<OsString>, env: I) -> Result<Prepassed, CliError>
where
    I: IntoIterator<Item = (OsString, OsString)>,
{
    let mut command = Args::command_with_error_table();
    command.build();
    let declared = declared_options(&command);
    let perl = perl_canonicals();
    let env = environment_map(env);

    let mut argv = argv.into_iter();
    let program = argv.next().unwrap_or_else(|| "vep".into());
    let options: Vec<OsString> = argv.collect();
    let tokens = spell(&options, &declared);
    if version_short_beside_arguments(&tokens) {
        return Err(CliError::cli(VERSION_SHORT_BESIDE_ARGUMENTS));
    }
    let asks_help_or_version = tokens
        .iter()
        .any(|t| matches!(t, Token::Option { long, .. } if long == "help" || long == "version"));
    let spelled: Vec<OsString> = tokens.iter().cloned().map(Token::into_os_string).collect();
    if asks_help_or_version {
        return Ok(Prepassed {
            argv: std::iter::once(program).chain(spelled).collect(),
            defined_off: Vec::new(),
        });
    }

    let named: BTreeSet<String> = tokens
        .iter()
        .filter_map(|t| match t {
            Token::Option { long, .. } => Some(long.clone()),
            _ => None,
        })
        .collect();
    let sources = sources_behind(
        command_line_value(&tokens, "config").as_deref(),
        command_line_value(&tokens, "dir").as_deref(),
        &named,
        &env,
        &declared,
        &perl,
    )?;
    let defined_off = defined_off_switches(&sources, &declared);
    let separator = spelled
        .iter()
        .position(|t| t == "--")
        .unwrap_or(spelled.len());
    let (own, rest) = spelled.split_at(separator);
    let mut argv = vec![program];
    argv.extend(lay_out(&sources, own, rest, &declared));
    Ok(Prepassed { argv, defined_off })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `@VEP_PARAMS`, Config.pm:73-276 at release/116.3, the comments dropped.
    const CONFIG_PM_PARAMS: &str = r#"our @VEP_PARAMS = (
  'help',
  'config=s',
  'input_file|i=s',
  'input_data|id=s',
  'format=s',
  'output_format=s',
  'delimiter=s',
  'no_check_variants_order',
  'species|s=s',
  'registry=s',
  'host=s',
  'port=s',
  'user|u=s',
  'password|pass=s',
  'db_version=i',
  'assembly|a=s',
  'grch37',
  'genomes',
  'refseq',
  'merged',
  'all_refseq',
  'gencode_basic',
  'gencode_primary',
  'is_multispecies=i',
  'transcript_filter=s@',
  'exclude_predicted',
  'minimal',
  'most_severe',
  'summary',
  'pick',
  'pick_allele',
  'per_gene',
  'pick_allele_gene',
  'flag_pick',
  'flag_pick_allele',
  'flag_pick_allele_gene',
  'flag_gencode_primary',
  'pick_order=s',
  'buffer_size=i',
  'failed=i',
  'gp',
  'chr=s',
  'check_ref',
  'lookup_ref',
  'check_existing',
  'check_svs',
  'no_check_alleles',
  'exclude_null_alleles',
  'check_frequency',
  'af',
  'af_1kg',
  'af_gnomade',
  'af_gnomad',
  'af_gnomadg',
  'old_maf',
  'max_af',
  'pubmed',
  'freq_filter=s',
  'freq_freq=f',
  'freq_gt_lt=s',
  'freq_pop=s',
  'filter_common',
  'allow_non_variant',
  'process_ref_homs',
  'individual=s',
  'individual_zyg=s',
  'phased',
  'fork=i',
  'dont_skip',
  'nearest=s',
  'distance=s',
  'clin_sig_allele=i',
  'clinvar_somatic_classification',
  'overlaps',
  'max_sv_size=i',
  'remove_hgvsp_version',
  'verbose|v',
  'quiet|q',
  'no_progress',
  'everything|e',
  'output_file|o=s@',
  'compress_output=s',
  'no_headers',
  'mask_header_cache_path',
  'stats_file|sf=s',
  'stats_text',
  'stats_html',
  'no_stats',
  'warning_file=s',
  'skipped_variants_file=s',
  'force_overwrite|force',
  'terms|t=s',
  'coding_only',
  'canonical',
  'mane',
  'mane_select',
  'tsl',
  'appris',
  'ccds',
  'xref_refseq',
  'uniprot',
  'protein',
  'protein_version',
  'biotype',
  'hgnc',
  'symbol',
  'transcript_version',
  'gene_version',
  'gene_phenotype',
  'mirna',
  'spdi',
  'ga4gh_vrs',
  'hgvs',
  'hgvsg',
  'hgvsg_use_accession',
  'hgvsp_use_prediction',
  'shift_hgvs=i',
  'ambiguous_hgvs',
  'sift=s',
  'polyphen=s',
  'humdiv',
  'condel=s',
  'variant_class',
  'regulatory',
  'cell_type=s',
  'regulatory_gff=s',
  'extended_promoters',
  'convert=s',
  'no_intergenic',
  'vcf',
  'solr',
  'json',
  'tab',
  'vcf_info_field=s',
  'keep_csq',
  'keep_ann',
  'lrg',
  'fields=s',
  'domains',
  'numbers',
  'total_length',
  'allele_number',
  'show_ref_allele',
  'no_escape',
  'ambiguity',
  'var_synonyms',
  'shift_3prime=i',
  'shift_genomic=i',
  'shift_length',
  'uploaded_allele',
  'database',
  'cache:s',
  'cache_version=i',
  'show_cache_info',
  'dir=s',
  'dir_cache=s',
  'full_cache_dir=s',
  'dir_plugins=s',
  'offline',
  'fasta|fa=s',
  'fasta_dir=s',
  'no_fasta',
  'sereal',
  'synonyms=s',
  'bam=s',
  'use_transcript_ref',
  'use_given_ref',
  'custom=s@',
  'tmpdir=s',
  'gff=s',
  'gtf=s',
  'bigwig=s',
  'phyloP=s@',
  'phastCons=s@',
  'ucsc_assembly=s',
  'ucsc_data_root=s',
  'custom_multi_allelic',
  'custom_suppress_filter',
  'plugin=s@',
  'safe',
  'debug',
);"#;

    fn command() -> Command {
        let mut command = Args::command_with_error_table();
        command.build();
        command
    }

    fn argv(args: &[&str]) -> Vec<OsString> {
        args.iter().map(OsString::from).collect()
    }

    fn strings(argv: &[OsString]) -> Vec<String> {
        argv.iter()
            .map(|a| a.to_string_lossy().into_owned())
            .collect()
    }

    fn env(pairs: &[(&str, &str)]) -> Vec<(OsString, OsString)> {
        pairs
            .iter()
            .map(|(k, v)| (OsString::from(k), OsString::from(v)))
            .collect()
    }

    /// A pre-pass with `home` as the home directory and `extra` in the
    /// environment, `vep` prepended to `args`.
    fn run(args: &[&str], home: &Path, extra: &[(&str, &str)]) -> Vec<String> {
        let mut pairs = vec![("HOME", home.to_str().unwrap())];
        pairs.extend(extra.iter().copied());
        let mut argv = vec!["vep"];
        argv.extend(args);
        strings(&prepass(self::argv(&argv), env(&pairs)).unwrap())
    }

    fn write(dir: &Path, name: &str, text: &str) -> String {
        let path = dir.join(name);
        fs::write(&path, text).unwrap();
        path.display().to_string()
    }

    fn owned(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    /// `VEP_PARAMS` is `@VEP_PARAMS`, entry for entry.
    #[test]
    fn vep_params_is_the_perl_list() {
        let perl: Vec<&str> = CONFIG_PM_PARAMS
            .lines()
            .filter_map(|l| l.trim().strip_prefix('\'')?.strip_suffix("',"))
            .collect();
        assert_eq!(perl.len(), 181);
        assert_eq!(VEP_PARAMS, perl.as_slice());
    }

    /// Getopt::Long's resolution over Ensembl VEP's names, each case checked
    /// against `GetOptions` with `@VEP_PARAMS`: exact in any case, a unique
    /// prefix, an alias, a single letter, an ambiguous or unknown prefix.
    #[test]
    fn resolve_reads_names_as_getopt_long_does() {
        let perl = perl_options();
        let name = |s: &str| resolve(s, &perl).map(|o| o.canonical.as_str());
        assert_eq!(name("format"), Some("format"));
        assert_eq!(name("FORMAT"), Some("format"));
        assert_eq!(name("form"), Some("format"));
        assert_eq!(name("fa"), Some("fasta"));
        assert_eq!(name("pass"), Some("password"));
        assert_eq!(name("i"), Some("input_file"));
        assert_eq!(name("js"), Some("json"));
        assert_eq!(name("af_gnomad"), Some("af_gnomad"));
        assert_eq!(name("no_s"), Some("no_stats"));
        assert_eq!(name("fo"), None);
        assert_eq!(name("hg"), None);
        assert_eq!(name("nonsense"), None);
        assert_eq!(name(""), None);
    }

    /// A `--config` file's flags come before the command line's, which win:
    /// the file's `format vcf` is dropped for the command line's `--format
    /// ensembl`, its `species`, `fork` and `offline` are added, its `quiet 0`
    /// leaves the switch off (`Config.pm:701`, `||=` at `:1101`).
    #[test]
    fn config_file_options_precede_command_line() {
        let dir = tempfile::tempdir().unwrap();
        let file = write(
            dir.path(),
            "flags.ini",
            "# a comment\nformat vcf\n\nspecies mus_musculus\nfork 4\noffline 1\nquiet 0\n",
        );
        let out = run(
            &["--config", &file, "--format", "ensembl", "-i", "x.vcf"],
            dir.path(),
            &[],
        );
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--fork=4",
                "--offline",
                "--species=mus_musculus",
                &format!("--config={file}"),
                "--format=ensembl",
                "--input_file=x.vcf",
            ])
        );
    }

    /// Keys of a configuration file as Ensembl VEP reads them, each case checked
    /// against `read_config_from_file` (`Config.pm:1077-1102`): a key with two
    /// dashes keeps one and names nothing; an unknown key, an alias, a prefix,
    /// another case and an indented line are dropped; a key with no value is
    /// dropped; a switch or `cache` set to `0` is off; a repeatable key takes
    /// the rest of its line and a double-quoted value keeps its blanks; a second
    /// value of a key is dropped unless the first was `0`; a key is kept as its
    /// `@VEP_PARAMS` spec takes it, a string option with its value (`cell_type
    /// HUVEC`); declared or not, a switch's `0` (`no_progress 0`) is off and a
    /// value's `0` (`clin_sig_allele 0`) is kept.
    #[test]
    fn config_unknown_key_is_dropped() {
        let dir = tempfile::tempdir().unwrap();
        let file = write(
            dir.path(),
            "flags.ini",
            "--dir_cache /double\n-dir_plugins /single\nno_such_flag 1\nfa ref.fa\ninput_f x.vcf\n\
             CACHE 1\n  format vcf\nno_headers\nquiet 0\nplugin \"Foo,bar baz\" x\n\
             plugin Second,x\ntranscript_filter biotype is protein_coding\n\
             species 'mus_musculus'\nbuffer_size 0\nbuffer_size 500\nfork 4\nfork 8\ncache 0\n\
             cell_type HUVEC\nphyloP x\nno_progress 0\nclin_sig_allele 0\n",
        );
        let out = run(&["--config", &file], dir.path(), &[]);
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--buffer_size=500",
                "--cell_type=HUVEC",
                "--clin_sig_allele=0",
                "--dir_plugins=/single",
                "--fork=4",
                "--species=mus_musculus",
                &format!("--config={file}"),
                "--phyloP=x",
                "--plugin=Foo,bar baz x",
                "--plugin=Second,x",
                "--transcript_filter=biotype is protein_coding",
            ])
        );
    }

    /// A source's `cache` as `Config.pm:716-718` reads it: `0` is off, as a
    /// switch's `0` is; the empty string, which Ensembl VEP's bare `--cache`
    /// also stores, is the flag alone; any other value is passed, and names the
    /// cache root when it is a directory. A required value is passed whatever
    /// its truth: `fork 0` is `--fork=0`.
    #[test]
    fn config_cache_value_is_read_as_perl_truth() {
        let dir = tempfile::tempdir().unwrap();
        for (line, added) in [
            ("cache 1\n", vec!["--cache=1"]),
            ("cache /data/vep\n", vec!["--cache=/data/vep"]),
            ("cache 0\n", vec![]),
            ("cache \"\"\n", vec!["--cache"]),
            ("fork 0\n", vec!["--fork=0"]),
        ] {
            let file = write(dir.path(), "flags.ini", line);
            let mut expected = owned(&["vep"]);
            expected.extend(owned(&added));
            expected.push(format!("--config={file}"));
            assert_eq!(
                run(&["--config", &file], dir.path(), &[]),
                expected,
                "{line:?}"
            );
        }
    }

    /// Getopt::Long takes an `=i` value matching `[-+]?_*[0-9][0-9_]*` and stores
    /// `0+$arg` with the underscores dropped; any other text is no integer.
    #[test]
    fn getopt_integer_reads_signs_zeros_and_underscores() {
        for (text, number) in [
            ("0", 0),
            ("1", 1),
            ("+1", 1),
            ("-1", -1),
            ("01", 1),
            ("00", 0),
            ("-0", 0),
            ("1_0", 10),
            ("_1", 1),
            ("1_", 1),
            ("-_1", -1),
            ("0_0", 0),
            ("1\n", 1),
            ("0\n", 0),
            ("9999999999999999999999", i64::MAX),
            ("-9999999999999999999999", -i64::MAX),
        ] {
            assert_eq!(getopt_integer(text), Some(number), "{text:?}");
        }
        for text in [
            "", "x", "1.5", "1e3", "_", "-", "+", "+-1", "_-1", " 1", "1 ", "0x1", "\n", "1\n\n",
            "\n1",
        ] {
            assert_eq!(getopt_integer(text), None, "{text:?}");
        }
    }

    /// A file's or the environment's `shift_3prime` and `shift_genomic` value
    /// stays the text its source holds (`Config.pm:1101, 1158`), true unless it
    /// is `0`, where the command line's is the number Getopt::Long stores: a
    /// text that is zero as a number is laid out as `1`, any other as written,
    /// and the command line's own value is left alone.
    #[test]
    fn a_lower_source_shift_value_is_read_as_text() {
        let dir = tempfile::tempdir().unwrap();
        for key in ["shift_3prime", "shift_genomic"] {
            for (value, added) in [
                ("00", "1"),
                ("+0", "1"),
                ("-0", "1"),
                ("0_0", "1"),
                ("0", "0"),
                ("2", "2"),
                ("01", "01"),
                ("x", "x"),
            ] {
                let token = format!("--{key}={added}");
                let file = write(dir.path(), "flags.ini", &format!("{key} {value}\n"));
                assert_eq!(
                    run(&["--config", &file], dir.path(), &[]),
                    owned(&["vep", &token, &format!("--config={file}")]),
                    "{key} {value}"
                );
                let variable = format!("VEP_{}", key.to_uppercase());
                assert_eq!(
                    run(&[], dir.path(), &[(&variable, value)]),
                    owned(&["vep", &token]),
                    "{variable}={value}"
                );
            }
            let flag = format!("--{key}=00");
            assert_eq!(run(&[&flag], dir.path(), &[]), owned(&["vep", &flag]));
        }
    }

    /// A lower source's `convert` value reaches the parser as written, `0` and
    /// the empty value included, so the deprecation check reads it as it reads
    /// the command line's (`Config.pm:977`). A value the command line gives,
    /// even a false one, keeps a lower source's from being read, where Ensembl
    /// VEP's `||=` replaces a false value by the lower source's true one
    /// (`Config.pm:1101, 1158`).
    #[test]
    fn a_lower_source_convert_value_is_laid_out_as_written() {
        let dir = tempfile::tempdir().unwrap();
        for (value, token) in [("0", "--convert=0"), ("vcf", "--convert=vcf")] {
            let file = write(dir.path(), "flags.ini", &format!("convert {value}\n"));
            assert_eq!(
                run(&["--config", &file], dir.path(), &[]),
                owned(&["vep", token, &format!("--config={file}")]),
                "convert {value}"
            );
            assert_eq!(
                run(&[], dir.path(), &[("VEP_CONVERT", value)]),
                owned(&["vep", token]),
                "VEP_CONVERT={value}"
            );
        }
        let file = write(dir.path(), "flags.ini", "convert \"\"\n");
        assert_eq!(
            run(&["--config", &file], dir.path(), &[]),
            owned(&["vep", "--convert=", &format!("--config={file}")])
        );
        let file = write(dir.path(), "flags.ini", "convert vcf\n");
        assert_eq!(
            run(
                &["--convert", "0", "--config", &file],
                dir.path(),
                &[("VEP_CONVERT", "vcf")]
            ),
            owned(&["vep", "--convert=0", &format!("--config={file}")])
        );
    }

    /// A key that is an Ensembl VEP option of its own but an alias here
    /// (`full_cache_dir`, `Config.pm:243`, the alias of `--json_cache`) names
    /// the option in a file and in the environment, as Ensembl VEP stores and
    /// reads it, and keeps that name on the token, the key Ensembl VEP prints in
    /// its command-line header (`Config.pm:734-745`), from every source and in
    /// any spelling of the command line's own; the long name stays its own; an
    /// alias both treat as one (`fa`) still names nothing.
    #[test]
    fn source_key_that_is_an_ensembl_vep_option_names_its_alias_target() {
        let dir = tempfile::tempdir().unwrap();
        let file = write(dir.path(), "flags.ini", "full_cache_dir /x\nfa ref.fa\n");
        assert_eq!(
            run(&["--config", &file], dir.path(), &[]),
            owned(&["vep", "--full_cache_dir=/x", &format!("--config={file}")])
        );
        assert_eq!(
            run(
                &["-i", "x.vcf"],
                dir.path(),
                &[("VEP_FULL_CACHE_DIR", "/y"), ("VEP_FA", "ref.fa")]
            ),
            owned(&["vep", "--full_cache_dir=/y", "--input_file=x.vcf"])
        );
        for spelling in [
            "--full_cache_dir",
            "-full_cache_dir",
            "--FULL_CACHE_DIR",
            "--full_cache",
        ] {
            assert_eq!(
                run(
                    &[spelling, "/z"],
                    dir.path(),
                    &[("VEP_FULL_CACHE_DIR", "/y")]
                ),
                owned(&["vep", "--full_cache_dir=/z"]),
                "{spelling}"
            );
        }
        assert_eq!(
            run(
                &["--json_cache", "/z"],
                dir.path(),
                &[("VEP_FULL_CACHE_DIR", "/y")]
            ),
            owned(&["vep", "--json_cache=/z"])
        );
        let file = write(dir.path(), "long.ini", "json_cache /x\n");
        assert_eq!(
            run(&["--config", &file], dir.path(), &[]),
            owned(&["vep", "--json_cache=/x", &format!("--config={file}")])
        );
    }

    /// A `--config` file that cannot be opened is fatal with Ensembl VEP's text
    /// (`Config.pm:1075`), the path as given.
    #[test]
    fn config_missing_file_is_fatal() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path().to_str().unwrap();
        let missing = dir.path().join("absent.ini").display().to_string();
        let err = prepass(
            argv(&["vep", "--config", &missing, "-i", "x.vcf"]),
            env(&[("HOME", home)]),
        )
        .unwrap_err();
        assert_eq!(
            err.to_string(),
            format!("ERROR: [config] Could not open config file \"{missing}\"")
        );
        let err = prepass(
            argv(&["vep", "--config=~/absent.ini"]),
            env(&[("HOME", home)]),
        )
        .unwrap_err();
        assert_eq!(
            err.to_string(),
            "ERROR: [config] Could not open config file \"~/absent.ini\""
        );
    }

    /// `<dir>/vep.ini` is read when it is a file (`Config.pm:704-705`): under
    /// `--dir`, under the `--config` file's `dir`, else under `~/.vep`; its
    /// flags come after the `--config` file's and before the command line's.
    #[test]
    fn vep_ini_under_dir_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path();
        fs::create_dir_all(home.join(".vep")).unwrap();
        write(&home.join(".vep"), INI_FILE, "symbol 1\nformat vcf\n");
        assert_eq!(
            run(&["-i", "x.vcf"], home, &[]),
            owned(&["vep", "--format=vcf", "--symbol", "--input_file=x.vcf"])
        );

        let other = home.join("other");
        fs::create_dir_all(&other).unwrap();
        write(&other, INI_FILE, "biotype 1\n");
        let other_dir = other.display().to_string();
        let dir_flag = format!("--dir={other_dir}");
        assert_eq!(
            run(&["--dir", &other_dir, "-i", "x.vcf"], home, &[]),
            owned(&["vep", "--biotype", &dir_flag, "--input_file=x.vcf"])
        );

        let config = write(
            home,
            "flags.ini",
            &format!("dir {other_dir}\nformat ensembl\n"),
        );
        assert_eq!(
            run(&["--config", &config, "-i", "x.vcf"], home, &[]),
            owned(&[
                "vep",
                "--biotype",
                &dir_flag,
                "--format=ensembl",
                &format!("--config={config}"),
                "--input_file=x.vcf",
            ])
        );

        let empty = home.join("empty");
        fs::create_dir_all(&empty).unwrap();
        let empty = empty.display().to_string();
        assert_eq!(
            run(&["--dir", &empty], home, &[]),
            owned(&["vep", &format!("--dir={empty}")])
        );

        // A directory named vep.ini reads as nothing (Perl's open on a directory
        // yields no line in Ensembl VEP), and a --dir that reads as false is the
        // default (`$config->{dir} || $DEFAULTS{dir}`).
        let as_dir = home.join("as_dir");
        fs::create_dir_all(as_dir.join(INI_FILE)).unwrap();
        let as_dir = as_dir.display().to_string();
        assert_eq!(
            run(&["--dir", &as_dir], home, &[]),
            owned(&["vep", &format!("--dir={as_dir}")])
        );
        assert_eq!(
            run(&["--dir", "", "-i", "x.vcf"], home, &[]),
            owned(&[
                "vep",
                "--format=vcf",
                "--symbol",
                "--dir=",
                "--input_file=x.vcf"
            ])
        );
    }

    /// A leading `~` is the home directory; without `HOME` the default `~/.vep`
    /// is `/.vep`, as `$ENV{HOME}.'/.vep'` is with the variable unset
    /// (`Config.pm:294`), so nothing is read relative to the working directory.
    #[test]
    fn tilde_expands_to_home_or_to_the_root() {
        let expand = |path: &str, home: &str| expand_tilde(OsStr::new(path), home);
        assert_eq!(expand("~/.vep", "/home/u"), PathBuf::from("/home/u/.vep"));
        assert_eq!(expand("~", "/home/u"), PathBuf::from("/home/u"));
        assert_eq!(expand(DEFAULT_DIR, ""), PathBuf::from("/.vep"));
        assert_eq!(expand("~", ""), PathBuf::from("/"));
        assert_eq!(expand("/abs/dir", "/home/u"), PathBuf::from("/abs/dir"));
        assert_eq!(expand("rel/~/dir", "/home/u"), PathBuf::from("rel/~/dir"));
    }

    /// A configuration file is read as bytes: a comment or a value holding a byte
    /// that is not UTF-8 does not make the file unreadable, as it does not for
    /// Ensembl VEP.
    #[test]
    fn config_file_bytes_outside_utf8_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("latin1.ini");
        fs::write(&path, b"# caf\xe9\noffline 1\nspecies caf\xe9\n").unwrap();
        let file = path.display().to_string();
        assert_eq!(
            run(&["--config", &file], dir.path(), &[]),
            owned(&[
                "vep",
                "--offline",
                "--species=caf\u{FFFD}",
                &format!("--config={file}"),
            ])
        );
    }

    /// The command line beats the `--config` file, which beats `vep.ini`, which
    /// beats the environment (`Config.pm:693-698`); a repeatable flag collects
    /// every source's values in that order (`:757`, `:1097`, `:1152`).
    #[test]
    fn priority_order_command_line_config_ini_env() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path();
        fs::create_dir_all(home.join(".vep")).unwrap();
        write(
            &home.join(".vep"),
            INI_FILE,
            "format ensembl\nspecies ini_species\nfork 2\nplugin FromIni\n",
        );
        let config = write(
            home,
            "flags.ini",
            "format vcf\nfork 4\nbuffer_size 100\nplugin FromConfig\n",
        );
        let out = run(
            &[
                &format!("--config={config}"),
                "--format=id",
                "--plugin=FromCommandLine",
            ],
            home,
            &[
                ("VEP_FORMAT", "spdi"),
                ("VEP_SPECIES", "env_species"),
                ("VEP_FORK", "8"),
                ("VEP_BUFFER_SIZE", "7"),
                ("VEP_DISTANCE", "9"),
                ("VEP_PLUGIN", "FromEnv"),
            ],
        );
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--distance=9",
                "--species=ini_species",
                "--buffer_size=100",
                "--fork=4",
                &format!("--config={config}"),
                "--format=id",
                "--plugin=FromCommandLine",
                "--plugin=FromConfig",
                "--plugin=FromIni",
                "--plugin=FromEnv",
            ])
        );
    }

    /// A `--config` path that is a directory reads as nothing, as Perl's
    /// two-argument `open` of a directory succeeds and its readline yields no
    /// line (`Config.pm:1075-1077`), the run going on; a path Perl reads as
    /// false, `0` or the empty value, is refused before any open with the
    /// empty string shown (`:1073`).
    #[test]
    fn config_directory_reads_nothing_and_a_false_path_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path();
        let as_dir = home.join("flags.d");
        fs::create_dir_all(&as_dir).unwrap();
        let as_dir = as_dir.display().to_string();
        assert_eq!(
            run(&["--config", &as_dir, "-i", "x.vcf"], home, &[]),
            owned(&["vep", &format!("--config={as_dir}"), "--input_file=x.vcf"])
        );
        let tilde_dir = home.join("tilde.d");
        fs::create_dir_all(&tilde_dir).unwrap();
        assert_eq!(
            run(&["--config=~/tilde.d"], home, &[]),
            owned(&["vep", "--config=~/tilde.d"])
        );

        for false_path in ["0", ""] {
            let err = prepass(
                argv(&["vep", "--config", false_path, "-i", "x.vcf"]),
                env(&[("HOME", home.to_str().unwrap())]),
            )
            .unwrap_err();
            assert_eq!(
                err.to_string(),
                "ERROR: [config] Could not open config file \"\"",
                "--config {false_path:?}"
            );
        }
    }

    /// The sources behind a companion tool's command line are the ones behind
    /// `vep`'s, read for `variant_recoder` by the same `Config->new`
    /// (`Config.pm:701-708`): the `--config` file, then `<dir>/vep.ini`, under
    /// `~/.vep` or under the file's `dir`, then the `VEP_*` variables; the
    /// command line beats the file, the file the ini, the ini the environment,
    /// and a repeatable flag collects every source's values in that order. A
    /// key the tool has no option for is dropped from every source, as is one
    /// of the options the tool presets or reads before the sources (`fields`,
    /// `buffer_size`, `ambiguous_hgvs`, `pretty`, `var_synonyms`:
    /// `VariantRecoder.pm:96-152`), a false switch adds nothing, and `--help`
    /// reads no source.
    #[test]
    fn companion_reads_vep_ini_and_the_environment_behind_the_config_file() {
        use crate::recode_cli::{RecodeArgs, PRESET_OPTIONS};
        use clap::CommandFactory;

        let dir = tempfile::tempdir().unwrap();
        let home = dir.path();
        fs::create_dir_all(home.join(".vep")).unwrap();
        write(
            &home.join(".vep"),
            INI_FILE,
            "input_data ini_input\nspecies ini_species\nshift_hgvs 0\nfields ini_fields\n\
             buffer_size 2\nambiguous_hgvs 1\ndir_plugins /ini/plugins\n\
             transcript_filter from_ini\n",
        );
        let config = write(
            home,
            "recode.ini",
            "input_data config_input\nshift_hgvs 1\nminimal 1\npretty 1\nconfig other.ini\n\
             transcript_filter from_config\n",
        );
        let mut command = RecodeArgs::command();
        command.build();
        let companion = |args: &[&str], extra: &[(&str, &str)]| {
            let mut pairs = vec![("HOME", home.to_str().unwrap())];
            pairs.extend(extra.iter().copied());
            let spelled = spell_as_declared(&argv(args), &command);
            strings(&with_config_file(&command, &spelled, env(&pairs), PRESET_OPTIONS).unwrap())
        };
        let variables = [
            ("VEP_INPUT_DATA", "env_input"),
            ("VEP_SPECIES", "env_species"),
            ("VEP_SHIFT_HGVS", "0"),
            ("VEP_ASSEMBLY", "GRCh38"),
            ("VEP_FIELDS", "env_fields"),
            ("VEP_VAR_SYNONYMS", "1"),
            ("VEP_DIR_PLUGINS", "/env/plugins"),
            ("VEP_MINIMAL", "0"),
            ("VEP_TRANSCRIPT_FILTER", "from_env"),
        ];
        assert_eq!(
            companion(
                &[
                    "--config",
                    &config,
                    "--id",
                    "cmd_input",
                    "--fields",
                    "hgvsg"
                ],
                &variables
            ),
            owned(&[
                "--assembly=GRCh38",
                "--species=ini_species",
                "--minimal",
                "--shift_hgvs=1",
                "--config",
                &config,
                "--input_data",
                "cmd_input",
                "--fields",
                "hgvsg",
                "--transcript_filter=from_config",
                "--transcript_filter=from_ini",
                "--transcript_filter=from_env",
            ])
        );
        // No `--config`: the ini and the environment alone.
        assert_eq!(
            companion(&["--fields", "hgvsg"], &variables),
            owned(&[
                "--assembly=GRCh38",
                "--input_data=ini_input",
                "--shift_hgvs=0",
                "--species=ini_species",
                "--fields",
                "hgvsg",
                "--transcript_filter=from_ini",
                "--transcript_filter=from_env",
            ])
        );
        // The file's `dir` places `vep.ini`, though the tool has no `--dir`.
        let other = home.join("other");
        fs::create_dir_all(&other).unwrap();
        write(&other, INI_FILE, "species other_species\n");
        let placed = write(home, "placed.ini", &format!("dir {}\n", other.display()));
        assert_eq!(
            companion(&["--config", &placed, "--id", "rs1"], &[]),
            owned(&[
                "--species=other_species",
                "--config",
                &placed,
                "--input_data",
                "rs1",
            ])
        );
        assert_eq!(
            companion(&["--help", "--config", "/no/such.ini"], &variables),
            owned(&["--help", "--config", "/no/such.ini"])
        );
    }

    /// `VEP_*` variables as `read_config_from_environment` reads them, each case
    /// checked against it (`Config.pm:1131-1161`): the name lower-cased, an
    /// empty value skipped, an unknown name, an alias and a prefix dropped,
    /// `VEP_FASTA_STORAGE` left to the FASTA reader, a switch set to `0` off, a
    /// variable without the upper-case `VEP_` prefix not read.
    #[test]
    fn environment_keys_follow_read_config_from_environment() {
        let dir = tempfile::tempdir().unwrap();
        let out = run(
            &["-i", "x.vcf"],
            dir.path(),
            &[
                ("VEP_Cache_Version", "115"),
                ("VEP_OFFLINE", "1"),
                ("VEP_QUIET", ""),
                ("VEP_NO_HEADERS", "0"),
                ("VEP_NO_SUCH_FLAG", "1"),
                ("VEP_FA", "ref.fa"),
                ("VEP_INPUT_F", "y.vcf"),
                ("VEP_FASTA_STORAGE", "memory"),
                ("VEP_CONFIG", "other.ini"),
                ("VEP_HELP", "1"),
                ("vep_symbol", "1"),
                ("VEPSYMBOL", "1"),
            ],
        );
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--cache_version=115",
                "--offline",
                "--input_file=x.vcf"
            ])
        );
    }

    /// No variable names a flag spelled with a capital letter: `VEP_PHYLOP`,
    /// `VEP_phyloP` and `VEP_PHASTCONS` are read as `phylop` and `phastcons`,
    /// which name nothing, as Ensembl VEP stores them under those keys
    /// (`Config.pm:1140, 1158`) where its option sets read `phyloP` and
    /// `phastCons` (`:553, 560`). `VEP_UCSC_ASSEMBLY` beside them sets its
    /// flag, and `VEP_CUSTOM` adds the source such a flag would name.
    #[test]
    fn environment_names_no_flag_spelled_with_a_capital_letter() {
        let dir = tempfile::tempdir().unwrap();
        let out = run(
            &["-i", "x.vcf"],
            dir.path(),
            &[
                ("VEP_PHYLOP", "100"),
                ("VEP_phyloP", "20"),
                ("VEP_PHASTCONS", "30"),
                ("VEP_UCSC_ASSEMBLY", "hg38"),
                ("VEP_CUSTOM", "scores.bw,phyloP100way,bigwig,exact"),
            ],
        );
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--ucsc_assembly=hg38",
                "--input_file=x.vcf",
                "--custom=scores.bw,phyloP100way,bigwig,exact",
            ])
        );
    }

    /// The command line's own spellings: one dash before a long name or a unique
    /// prefix (`-form`, `-sa`), any case, `=` values, the value after a flag left
    /// alone, a short cluster and a name no option has left as written, everything
    /// after `--` left as written. A name the parser's flags lack or leave
    /// ambiguous is written as Getopt::Long reads it over `@VEP_PARAMS`: an
    /// option the parser declares by its long name (`-cell_type` and `-ve` are
    /// `cell_type` and `verbose`; `--js` and `--pa` are `json` and `password`,
    /// which vep-rs's own flags make ambiguous; `-id` is the `input_data`
    /// alias). The error table's `-af_exac` is written in full for
    /// its refusal. A short with its value attached is the long form; a one-dash
    /// prefix of several options (`-sh`, `-ou`) is written with two dashes for
    /// clap to refuse, never as a short with a value.
    #[test]
    fn command_line_spellings_are_normalised() {
        let dir = tempfile::tempdir().unwrap();
        let out = run(
            &[
                "-cache",
                "-offline",
                "-form",
                "vcf",
                "--FORMAT=vcf",
                "-fa",
                "ref.fa",
                "--input_file=-cache",
                "-i",
                "-offline",
                "--max_sv_size=-1",
                "-qv",
                "--fo",
                "-nonsense",
                "-cell_type",
                "HUVEC",
                "-sa",
                "-af_exac",
                "-ve",
                "--js",
                "--pa",
                "secret",
                "-ob.txt",
                "-o=b=c.txt",
                "-id",
                "rs1",
                "-sh",
                "-ou",
                "--",
                "-cache",
            ],
            dir.path(),
            &[],
        );
        assert_eq!(
            out,
            owned(&[
                "vep",
                "--cache",
                "--offline",
                "--format=vcf",
                "--format=vcf",
                "--fasta=ref.fa",
                "--input_file=-cache",
                "--input_file=-offline",
                "--max_sv_size=-1",
                "-qv",
                "--fo",
                "-nonsense",
                "--cell_type=HUVEC",
                "--safe",
                "--af_exac",
                "--verbose",
                "--json",
                "--password=secret",
                "--output_file=b.txt",
                "--output_file=b=c.txt",
                "--input_data=rs1",
                "--sh",
                "--ou",
                "--",
                "-cache",
            ])
        );
    }

    /// A source's repeatable flags come after the command line's own and before
    /// its `--` separator.
    #[test]
    fn source_lists_precede_the_separator() {
        let dir = tempfile::tempdir().unwrap();
        let out = run(
            &["--plugin", "A", "--", "x"],
            dir.path(),
            &[("VEP_PLUGIN", "B")],
        );
        assert_eq!(out, owned(&["vep", "--plugin=A", "--plugin=B", "--", "x"]));
    }

    /// A command line asking for the help or the version reads no source (vep:29
    /// exits before `Config.pm` opens a file), in any spelling.
    #[test]
    fn help_and_version_read_no_source() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path().to_str().unwrap();
        let missing = dir.path().join("absent.ini").display().to_string();
        let config = format!("--config={missing}");
        let cases: [(&[&str], &[&str]); 6] = [
            (&["--config", &missing, "-help"], &[&config, "--help"]),
            (&["--HELP", "--config", &missing], &["--help", &config]),
            (&["-h"], &["--help"]),
            (&["--vers"], &["--version"]),
            (
                &["--version", "--config", &missing],
                &["--version", &config],
            ),
            (&["-V"], &["--version"]),
        ];
        for (args, expected) in cases {
            let mut argv = vec!["vep"];
            argv.extend(args.iter().copied());
            let out = prepass(
                self::argv(&argv),
                env(&[("HOME", home), ("VEP_SYMBOL", "1")]),
            );
            let out = match out {
                Ok(out) => strings(&out),
                Err(e) => panic!("{args:?}: {e}"),
            };
            let mut want = vec!["vep"];
            want.extend(expected.iter().copied());
            assert_eq!(out, owned(&want), "{args:?}");
        }
    }

    /// The token after an option Ensembl VEP declares as a number (`fork=i`,
    /// `freq_freq=f`) is left unread when it asks for the help or the version, as
    /// Getopt::Long pushes back a value that is not a number, so `--fork --help`
    /// asks for the help and reads no source; so is the token after an option
    /// Ensembl VEP lacks (`pipeline_depth`, `writer_strategy`), which
    /// `pass_through` leaves in `@ARGV` with its would-be value unread; a string
    /// option takes the token as its value, as `=s` takes anything (`-i --help`
    /// names the file `--help`), an option of `@VEP_PARAMS` the parser lacks
    /// included (`bigwig=s`, Config.pm:262), and a number is a number's value.
    #[test]
    fn a_numeric_option_leaves_a_usage_token_unread() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path().to_str().unwrap();
        let cases: [(&[&str], &[&str]); 9] = [
            (&["--fork", "--help"], &["--fork", "--help"]),
            (&["--fork", "x", "--help"], &["--fork=x", "--help"]),
            (&["--freq_freq", "--version"], &["--freq_freq", "--version"]),
            (&["--buffer_size", "-help"], &["--buffer_size", "--help"]),
            (&["--fork", "4", "--help"], &["--fork=4", "--help"]),
            (
                &["--pipeline_depth", "--help"],
                &["--pipeline_depth", "--help"],
            ),
            (
                &["--writer_strategy", "--version"],
                &["--writer_strategy", "--version"],
            ),
            (
                &["--input_file", "--help"],
                &["--symbol", "--input_file=--help"],
            ),
            (&["--bigwig", "--help"], &["--symbol", "--bigwig=--help"]),
        ];
        for (args, expected) in cases {
            let mut argv = vec!["vep"];
            argv.extend(args.iter().copied());
            let out = prepass(
                self::argv(&argv),
                env(&[("HOME", home), ("VEP_SYMBOL", "1")]),
            )
            .unwrap_or_else(|e| panic!("{args:?}: {e}"));
            let mut want = vec!["vep"];
            want.extend(expected.iter().copied());
            assert_eq!(strings(&out), owned(&want), "{args:?}");
        }
    }

    /// `-V` is the version when it is the only argument. Beside any other token
    /// the line is refused before a source is read, since Getopt::Long's
    /// `ignore_case` reads the token as `-v` (`verbose|v`, Config.pm:159) and
    /// Ensembl VEP runs: a second `-V`, a `--help=1`, which asks for no help, and
    /// a number option, which leaves the token unread, are such tokens. A `--help`
    /// or `--version` spelled out wins, and a `-V` that is a string option's
    /// value, a letter of a cluster or behind `--` is not the version's short.
    #[test]
    fn version_short_beside_another_token_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let home = dir.path().to_str().unwrap();
        let missing = dir.path().join("absent.ini").display().to_string();
        let read = |args: &[&str]| {
            let mut argv = vec!["vep"];
            argv.extend(args.iter().copied());
            prepass(self::argv(&argv), env(&[("HOME", home)]))
                .map(|out| strings(&out))
                .map_err(|e| e.to_string())
        };
        let refused: [&[&str]; 9] = [
            &["-V", "-i", "in.vcf"],
            &["-i", "in.vcf", "-o", "out.txt", "-V"],
            &["-V", "--quiet"],
            &["-V", "-V"],
            &["--V", "--cache"],
            &["--fork", "-V"],
            &["-V", "--config", &missing],
            &["-V", "--help=1"],
            &["-V", "--", "x"],
        ];
        for args in refused {
            assert_eq!(
                read(args),
                Err(format!("ERROR: [cli] {VERSION_SHORT_BESIDE_ARGUMENTS}")),
                "{args:?}"
            );
        }
        let passed: [(&[&str], &[&str]); 7] = [
            (&["-V"], &["--version"]),
            (&["-V", "--help"], &["--version", "--help"]),
            (&["-h", "-V"], &["--help", "--version"]),
            (&["-V", "--version"], &["--version", "--version"]),
            (&["-i", "-V"], &["--input_file=-V"]),
            (&["-qV"], &["-qV"]),
            (&["--", "-V"], &["--", "-V"]),
        ];
        for (args, expected) in passed {
            let mut want = vec!["vep"];
            want.extend(expected.iter().copied());
            assert_eq!(read(args), Ok(owned(&want)), "{args:?}");
        }
    }

    /// The parser's options carry every spelling the pre-pass resolves: the long
    /// name, the aliases and the short; a short is matched by case, as clap
    /// reads it, so `-V` stays the version and `-v` the verbosity.
    #[test]
    fn parser_options_name_aliases_and_shorts() {
        let parser = parser_options(&command());
        let input = parser.iter().find(|o| o.canonical == "input_file").unwrap();
        assert_eq!(input.names, ["input_file", "i"]);
        let force = parser
            .iter()
            .find(|o| o.canonical == "force_overwrite")
            .unwrap();
        assert_eq!(force.names, ["force_overwrite", "force"]);
        let name = |s: &str| resolve(s, &parser).map(|o| o.canonical.as_str());
        assert_eq!(name("i"), Some("input_file"));
        assert_eq!(name("fa"), Some("fasta"));
        assert_eq!(name("full_cache_dir"), Some("json_cache"));
        assert_eq!(name("threads"), Some("fork"));
        assert_eq!(name("v"), Some("verbose"));
        assert_eq!(name("V"), Some("version"));
        assert_eq!(name("js"), None);
        let perl = perl_options();
        let known = known_options(&declared_options(&command()));
        let spelled = |s: &str| resolve_spelling(s, &known, &perl).map(|o| o.canonical.as_str());
        assert_eq!(spelled("js"), Some("json"));
        assert_eq!(spelled("json_c"), Some("json_cache"));
        assert_eq!(spelled("ver"), Some("verbose"));
        assert_eq!(spelled("vers"), Some("version"));
        assert_eq!(spelled("af_exac"), Some("af_exac"));
        assert_eq!(spelled("stats_h"), Some("stats_html"));
        assert_eq!(spelled("cell_t"), Some("cell_type"));
        assert_eq!(spelled("fo"), None);
        assert_eq!(spelled("qv"), None);
    }

    /// The directory `vep.ini` is read from by default is the parser's default
    /// for `--dir`, so the two cannot drift apart.
    #[test]
    fn default_dir_is_the_parsers_default() {
        let command = command();
        let dir = command
            .get_arguments()
            .find(|a| a.get_long() == Some("dir"))
            .unwrap();
        let defaults: Vec<String> = dir
            .get_default_values()
            .iter()
            .map(|v| v.to_string_lossy().into_owned())
            .collect();
        assert_eq!(defaults, [DEFAULT_DIR]);
    }

    /// Configuration lines split as `Config.pm:1082-1092` splits them, each case
    /// checked against that code: a run of blanks inside double quotes becomes
    /// one blank of the word, single quotes protect nothing, a leading blank
    /// makes an empty first word, trailing blanks nothing, and the first word
    /// keeps its quotes.
    #[test]
    fn config_lines_split_as_perl_does() {
        assert_eq!(split_config_line("format vcf"), ["format", "vcf"]);
        assert_eq!(
            split_config_line("plugin \"Foo bar, baz\" x  \"a  b\""),
            ["plugin", "Foo bar, baz", "x", "a b"]
        );
        assert_eq!(
            split_config_line("species 'mus musculus'"),
            ["species", "mus", "musculus"]
        );
        assert_eq!(split_config_line("  format vcf"), ["", "format", "vcf"]);
        assert_eq!(split_config_line("fork  4  "), ["fork", "4"]);
        assert_eq!(split_config_line("no_headers"), ["no_headers"]);
        assert_eq!(split_config_line("\"format\" vcf"), ["\"format\"", "vcf"]);
        assert_eq!(
            split_config_line("plugin \"Foo bar"),
            ["plugin \"Foo", "bar"]
        );
    }

    /// A companion tool's name that is a prefix of several options is written
    /// with two dashes whatever its own, with its attached value, so clap
    /// refuses the name instead of reading its letters as short forms. A
    /// prefix of one option resolves; a single letter, a name no option has,
    /// the value of an option and a token behind `--` stay as typed.
    #[test]
    fn a_companion_prefix_of_several_options_is_written_with_two_dashes() {
        let switch =
            |name: &'static str| clap::Arg::new(name).long(name).action(ArgAction::SetTrue);
        let valued = |name: &'static str| clap::Arg::new(name).long(name).num_args(1);
        let tool = Command::new("tool")
            .arg(switch("list").short('l'))
            .arg(valued("limit"))
            .arg(valued("input_file").short('i'))
            .arg(valued("pass"))
            .arg(valued("port"));
        let spelled = |args: &[&str]| strings(&spell_as_declared(&argv(args), &tool));
        assert_eq!(spelled(&["-li", "x"]), ["--li", "x"]);
        assert_eq!(spelled(&["--li", "x"]), ["--li", "x"]);
        assert_eq!(spelled(&["-li=3"]), ["--li=3"]);
        assert_eq!(spelled(&["-lis"]), ["--list"]);
        assert_eq!(spelled(&["-lim", "-li"]), ["--limit", "-li"]);
        assert_eq!(
            spelled(&["-p", "-lx", "--", "-li"]),
            ["-p", "-lx", "--", "-li"]
        );
    }
}
