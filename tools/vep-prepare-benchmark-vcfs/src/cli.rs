// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as CPython's `argparse` reads a parser of these nine options with
//! its defaults (`allow_abbrev`, no positionals): an option is named whole or by a prefix
//! that fits it alone, takes its value from the next argument or after `=`, and a later
//! value replaces an earlier one. `--canonical-contigs` takes the next argument as its value
//! only where that argument is not itself an option. Errors come in argparse's order: what
//! an option meets as it is read, then the required options that are missing, then the
//! arguments no option took.

use crate::pyfmt::{is_py_decimal, py_int, py_repr, py_stderr_os_str, PyInt};
use std::ffi::OsString;
use std::os::unix::ffi::OsStrExt;

pub const PROG: &str = "vep-prepare-benchmark-vcfs";

pub const USAGE: &str = "\
usage: vep-prepare-benchmark-vcfs [-h] --input-dir INPUT_DIR
                                  --output-dir OUTPUT_DIR
                                  [--chrom-mode {auto,keep,strip,add}]
                                  [--cache-dir CACHE_DIR]
                                  [--max-variants MAX_VARIANTS]
                                  [--assembly ASSEMBLY]
                                  [--canonical-contigs [FILE_OR_LIST]]
                                  [--no-canonical-contigs]";

pub const HELP: &str = "\
Prepare benchmark VCF inputs

options:
  -h, --help            show this help message and exit
  --input-dir INPUT_DIR
                        Directory containing benchmark .vcf.gz or .vcf files
  --output-dir OUTPUT_DIR
                        Directory to write prepared .vcf files
  --chrom-mode {auto,keep,strip,add}
                        Chromosome normalization mode (default: auto).
                        Independent of the canonical-contigs filter; if
                        --canonical-contigs is active, UCSC-style inputs are
                        auto-stripped regardless of this mode.
  --cache-dir CACHE_DIR
                        Converted JSON cache directory (used only by --chrom-
                        mode auto)
  --max-variants MAX_VARIANTS
                        If set, keep only first N variant lines per input file
  --assembly ASSEMBLY   Reference assembly metadata (recorded in manifest).
                        Absent, GRCh37 or GRCh38 selects the human canonical
                        set for a bare --canonical-contigs; another assembly's
                        set is named with --canonical-contigs <file or list>,
                        and a bare flag under it switches the filter off with
                        a warning.
  --canonical-contigs [FILE_OR_LIST]
                        Filter input VCFs to canonical contigs (DEFAULT).
                        Bare, the universal Ensembl-style human set
                        (1-22,X,Y,MT); with a value, the comma-separated names
                        given or the names in the file named, one per line or
                        comma-separated, `#` lines ignored (a cache's
                        info.json contig_order). UCSC-style chr-prefix inputs
                        are auto-stripped at filter time. Drops variants on
                        alt contigs (NT_*, NW_*, decoys), so that asymmetric
                        alt-contig handling between engines cannot inflate F1
                        divergences.
  --no-canonical-contigs
                        Opt out of canonical-contigs filtering and annotate
                        every contig the input carries, alt contigs included.";

/// The `--chrom-mode` a command line asks for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChromMode {
    Auto,
    Keep,
    Strip,
    Add,
}

impl ChromMode {
    const CHOICES: [(&'static str, ChromMode); 4] = [
        ("auto", ChromMode::Auto),
        ("keep", ChromMode::Keep),
        ("strip", ChromMode::Strip),
        ("add", ChromMode::Add),
    ];

    pub fn name(self) -> &'static str {
        match self {
            ChromMode::Auto => "auto",
            ChromMode::Keep => "keep",
            ChromMode::Strip => "strip",
            ChromMode::Add => "add",
        }
    }
}

/// What the two contig options leave: the bare flag (the default), a named set, or the
/// filter switched off.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ContigsFlag {
    Bare,
    Named(Vec<u8>),
    Off,
}

/// The options of one run, each value in the bytes the command line gave.
#[derive(Debug)]
pub struct Args {
    pub input_dir: Vec<u8>,
    pub output_dir: Vec<u8>,
    pub chrom_mode: ChromMode,
    pub cache_dir: Option<Vec<u8>>,
    pub max_variants: Option<PyInt>,
    pub assembly: Option<Vec<u8>>,
    pub canonical_contigs: ContigsFlag,
}

/// What a command line asks for.
#[derive(Debug)]
pub enum Command {
    Run(Args),
    Help,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Opt {
    Help,
    InputDir,
    OutputDir,
    ChromMode,
    CacheDir,
    MaxVariants,
    Assembly,
    CanonicalContigs,
    NoCanonicalContigs,
}

/// How many values an option takes.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Arity {
    None,
    One,
    Optional,
}

impl Opt {
    fn arity(self) -> Arity {
        match self {
            Opt::Help | Opt::NoCanonicalContigs => Arity::None,
            Opt::CanonicalContigs => Arity::Optional,
            _ => Arity::One,
        }
    }

    /// The name an error calls the option by.
    fn name(self) -> &'static str {
        match self {
            Opt::Help => "-h/--help",
            Opt::InputDir => "--input-dir",
            Opt::OutputDir => "--output-dir",
            Opt::ChromMode => "--chrom-mode",
            Opt::CacheDir => "--cache-dir",
            Opt::MaxVariants => "--max-variants",
            Opt::Assembly => "--assembly",
            Opt::CanonicalContigs => "--canonical-contigs",
            Opt::NoCanonicalContigs => "--no-canonical-contigs",
        }
    }
}

/// The option strings in the order a usage error lists them.
const OPTION_STRINGS: [(&str, Opt); 10] = [
    ("-h", Opt::Help),
    ("--help", Opt::Help),
    ("--input-dir", Opt::InputDir),
    ("--output-dir", Opt::OutputDir),
    ("--chrom-mode", Opt::ChromMode),
    ("--cache-dir", Opt::CacheDir),
    ("--max-variants", Opt::MaxVariants),
    ("--assembly", Opt::Assembly),
    ("--canonical-contigs", Opt::CanonicalContigs),
    ("--no-canonical-contigs", Opt::NoCanonicalContigs),
];

/// One way to read an argument as an option: the option string, the text the argument
/// carries for it, and whether an `=` introduced that text.
#[derive(Clone, Copy)]
struct Candidate<'a> {
    string: &'static str,
    opt: Opt,
    explicit: Option<&'a [u8]>,
    after_equals: bool,
}

/// How one argument reads before any is consumed.
enum Reading<'a> {
    /// A value for the option before it, or a stray positional.
    Argument,
    /// The options the argument can name; more than one is an ambiguous prefix.
    Options(Vec<Candidate<'a>>),
    /// Shaped like an option and naming none.
    Unknown,
    /// The `--` marker, or an argument after it.
    Marker,
}

fn exact(string: &[u8]) -> Option<(&'static str, Opt)> {
    OPTION_STRINGS
        .iter()
        .copied()
        .find(|(option, _)| option.as_bytes() == string)
}

fn candidate<'a>(
    (string, opt): (&'static str, Opt),
    explicit: Option<&'a [u8]>,
    after_equals: bool,
) -> Candidate<'a> {
    Candidate {
        string,
        opt,
        explicit,
        after_equals,
    }
}

/// Reads one argument that stands before any `--` marker. An argument that starts with `-`
/// is an option unless it is `-` alone, starts like a negative number (`-` or `-.` and a
/// decimal digit of any script) or holds a space; such a one is a value.
fn read(arg: &[u8]) -> Reading<'_> {
    if arg.first() != Some(&b'-') {
        return Reading::Argument;
    }
    if let Some(found) = exact(arg) {
        return Reading::Options(vec![candidate(found, None, false)]);
    }
    if arg.len() == 1 {
        return Reading::Argument;
    }
    let (head, inline) = match memchr::memchr(b'=', arg) {
        Some(i) => (&arg[..i], Some(&arg[i + 1..])),
        None => (arg, None),
    };
    if let (Some(found), Some(_)) = (exact(head), inline) {
        return Reading::Options(vec![candidate(found, inline, true)]);
    }
    let mut candidates = Vec::new();
    for found in OPTION_STRINGS {
        let string = found.0.as_bytes();
        if arg[1] == b'-' {
            if string.starts_with(head) {
                candidates.push(candidate(found, inline, inline.is_some()));
            }
        } else if string == &arg[..2] {
            // A one-letter option with more characters behind it.
            candidates.push(candidate(found, Some(&arg[2..]), false));
        } else if string.starts_with(head) {
            candidates.push(candidate(found, inline, inline.is_some()));
        }
    }
    if !candidates.is_empty() {
        return Reading::Options(candidates);
    }
    let after_dash = arg[1..].strip_prefix(b".").unwrap_or(&arg[1..]);
    let negative_number = (1..=after_dash.len().min(4))
        .find_map(|length| std::str::from_utf8(&after_dash[..length]).ok())
        .and_then(|text| text.chars().next())
        .is_some_and(is_py_decimal);
    if negative_number || arg.contains(&b' ') {
        Reading::Argument
    } else {
        Reading::Unknown
    }
}

fn shown(bytes: &[u8]) -> String {
    py_stderr_os_str(bytes)
}

fn ignored_explicit(opt: Opt, value: &[u8]) -> String {
    format!(
        "argument {}: ignored explicit argument {}",
        opt.name(),
        py_repr(value)
    )
}

/// The values the options have set so far.
struct Values {
    input_dir: Option<Vec<u8>>,
    output_dir: Option<Vec<u8>>,
    chrom_mode: ChromMode,
    cache_dir: Option<Vec<u8>>,
    max_variants: Option<PyInt>,
    assembly: Option<Vec<u8>>,
    canonical_contigs: ContigsFlag,
}

impl Values {
    /// Sets what one option carries. `Ok(true)` is a request for the help text.
    fn take(&mut self, opt: Opt, value: Option<&[u8]>) -> Result<bool, String> {
        let owned = || value.map(<[u8]>::to_vec);
        match opt {
            Opt::Help => return Ok(true),
            Opt::InputDir => self.input_dir = owned(),
            Opt::OutputDir => self.output_dir = owned(),
            Opt::CacheDir => self.cache_dir = owned(),
            Opt::Assembly => self.assembly = owned(),
            Opt::ChromMode => {
                let value = value.unwrap_or_default();
                self.chrom_mode = ChromMode::CHOICES
                    .iter()
                    .find(|(name, _)| name.as_bytes() == value)
                    .map(|(_, mode)| *mode)
                    .ok_or_else(|| {
                        format!(
                            "argument --chrom-mode: invalid choice: {} (choose from 'auto', 'keep', 'strip', 'add')",
                            py_repr(value)
                        )
                    })?;
            }
            Opt::MaxVariants => {
                let value = value.unwrap_or_default();
                let read = std::str::from_utf8(value).ok().and_then(py_int);
                self.max_variants = Some(read.ok_or_else(|| {
                    format!(
                        "argument --max-variants: invalid int value: {}",
                        py_repr(value)
                    )
                })?);
            }
            Opt::CanonicalContigs => {
                self.canonical_contigs = match value {
                    Some(value) => ContigsFlag::Named(value.to_vec()),
                    None => ContigsFlag::Bare,
                };
            }
            Opt::NoCanonicalContigs => self.canonical_contigs = ContigsFlag::Off,
        }
        Ok(false)
    }
}

pub fn parse(argv: &[OsString]) -> Result<Command, String> {
    let args: Vec<&[u8]> = argv.iter().map(|arg| arg.as_bytes()).collect();
    let marker = args.iter().position(|arg| *arg == b"--");
    let readings: Vec<Reading> = args
        .iter()
        .enumerate()
        .map(|(index, arg)| match marker {
            Some(at) if index >= at => Reading::Marker,
            _ => read(arg),
        })
        .collect();
    let mut values = Values {
        input_dir: None,
        output_dir: None,
        chrom_mode: ChromMode::Auto,
        cache_dir: None,
        max_variants: None,
        assembly: None,
        canonical_contigs: ContigsFlag::Bare,
    };
    let mut unread: Vec<String> = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let arg = args[index];
        let candidates = match &readings[index] {
            Reading::Options(candidates) => candidates,
            Reading::Argument | Reading::Unknown | Reading::Marker => {
                unread.push(shown(arg));
                index += 1;
                continue;
            }
        };
        index += 1;
        if candidates.len() > 1 {
            let matches: Vec<&str> = candidates.iter().map(|c| c.string).collect();
            return Err(format!(
                "ambiguous option: {} could match {}",
                shown(arg),
                matches.join(", ")
            ));
        }
        let Candidate {
            string,
            mut opt,
            mut explicit,
            mut after_equals,
        } = candidates[0];
        let one_letter = string.len() == 2;
        // What the argument sets, in order: more than one where one-letter options run
        // together behind a single `-`.
        let mut taken: Vec<(Opt, Option<&[u8]>)> = Vec::new();
        loop {
            let Some(value) = explicit else {
                let next_is_value = matches!(readings.get(index), Some(Reading::Argument));
                match opt.arity() {
                    Arity::None => taken.push((opt, None)),
                    Arity::One if next_is_value => {
                        taken.push((opt, Some(args[index])));
                        index += 1;
                    }
                    Arity::One => {
                        return Err(format!("argument {}: expected one argument", opt.name()))
                    }
                    Arity::Optional if next_is_value => {
                        taken.push((opt, Some(args[index])));
                        index += 1;
                    }
                    Arity::Optional => taken.push((opt, None)),
                }
                break;
            };
            if opt.arity() != Arity::None {
                taken.push((opt, Some(value)));
                break;
            }
            if !one_letter || value.is_empty() || after_equals || value[0] == b'-' {
                return Err(ignored_explicit(opt, value));
            }
            // The characters behind a one-letter option that takes no value read as more
            // one-letter options; one that names no option is left unread.
            taken.push((opt, None));
            if value[0] != b'h' {
                unread.push(format!("-{}", shown(value)));
                break;
            }
            opt = Opt::Help;
            (explicit, after_equals) = match &value[1..] {
                [] => (None, false),
                [b'=', rest @ ..] => (Some(rest), true),
                rest => (Some(rest), false),
            };
        }
        for (opt, value) in taken {
            if values.take(opt, value)? {
                return Ok(Command::Help);
            }
        }
    }
    let missing: Vec<&str> = [
        ("--input-dir", values.input_dir.is_none()),
        ("--output-dir", values.output_dir.is_none()),
    ]
    .iter()
    .filter(|(_, absent)| *absent)
    .map(|(name, _)| *name)
    .collect();
    match (values.input_dir, values.output_dir) {
        (Some(input_dir), Some(output_dir)) if unread.is_empty() => Ok(Command::Run(Args {
            input_dir,
            output_dir,
            chrom_mode: values.chrom_mode,
            cache_dir: values.cache_dir,
            max_variants: values.max_variants,
            assembly: values.assembly,
            canonical_contigs: values.canonical_contigs,
        })),
        (Some(_), Some(_)) => Err(format!("unrecognized arguments: {}", unread.join(" "))),
        _ => Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )),
    }
}

#[cfg(test)]
#[path = "cli_parse_tests.rs"]
mod parse_tests;
