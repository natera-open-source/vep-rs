// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as CPython's `argparse` reads a parser of these options with its
//! defaults (`allow_abbrev`, no positionals): an option is named whole or by a prefix that
//! fits it alone, takes its value from the next argument or after `=`, and a later value
//! replaces an earlier one unless the option is repeatable.

use crate::pyfmt::py_int_clamped;
use crate::pypath;

pub const PROGRAM: &str = "vep-divergence-reevaluate";
pub const USAGE: &str = "\
usage: vep-divergence-reevaluate [-h] [--image IMAGE]
                                 [--ensembl-vep ENSEMBL_VEP] [--vep VEP]
                                 [--perl-cache-dir PERL_CACHE_DIR]
                                 [--plugins-dir PLUGINS_DIR]
                                 [--tool-dir TOOL_DIRS]
                                 [--extra-dir NAME=PATH] [--entry ENTRY]
                                 [--jobs JOBS] [--work-dir WORK_DIR] [--keep]
                                 [--timeout TIMEOUT] [--expect-all-present]
                                 [--json JSON] [--check]
                                 [--reference REFERENCE] [--page PAGE]
                                 [--reproductions REPRODUCTIONS] [--list]";
pub const HELP: &str = "\
Re-run every intended divergence on Ensembl VEP and vep-rs and say whether
each is still present.

options:
  -h, --help            show this help message and exit
  --image IMAGE         the Ensembl VEP image, a tag or a digest
  --ensembl-vep ENSEMBL_VEP
                        a local Ensembl VEP command to run in place of the
                        image (`vep`, `filter_vep` and `variant_recoder`
                        beside it)
  --vep VEP             the vep-rs binary
  --perl-cache-dir PERL_CACHE_DIR
                        a root holding
                        <species>[_merged|_refseq]/<version>_<assembly>/
                        (repeatable, searched in order)
  --plugins-dir PLUGINS_DIR
                        the Ensembl VEP plugin modules (the VEP_plugins
                        checkout)
  --tool-dir TOOL_DIRS  a directory holding a tool a reproduction's run.tools
                        names (repeatable); put first on Ensembl VEP's PATH
  --extra-dir NAME=PATH
                        a directory a reproduction names as {NAME}
                        (repeatable), such as a plugin repository outside the
                        plugin modules
  --entry ENTRY         run this entry alone (repeatable)
  --jobs JOBS           entries run at once
  --work-dir WORK_DIR   where the engines write (default: a temporary
                        directory)
  --keep                keep the engines' outputs
  --timeout TIMEOUT     seconds allowed for one engine run
  --expect-all-present  exit 1 unless every entry is still-present
  --json JSON           write the results as JSON here
  --check               check the page against the reproductions and run
                        nothing
  --reference REFERENCE
                        the Ensembl VEP release every entry is verified
                        against (--check)
  --page PAGE           the index page (--check)
  --reproductions REPRODUCTIONS
                        the reproduction directory
  --list                list the reproductions and run nothing";

#[derive(Clone, Debug, PartialEq)]
pub struct Args {
    pub image: Option<String>,
    pub ensembl_vep: Option<String>,
    pub vep: Option<String>,
    /// Whether the `--vep` value was written with a directory (`./vep`, `bin/vep`). The
    /// value itself does not show it for `./vep`, which is kept as `vep`.
    pub vep_names_a_directory: bool,
    pub perl_cache_dir: Vec<String>,
    pub plugins_dir: Option<String>,
    pub tool_dirs: Vec<String>,
    pub extra_dir: Vec<String>,
    pub entry: Vec<String>,
    pub jobs: i64,
    pub work_dir: Option<String>,
    pub keep: bool,
    pub timeout: i64,
    pub expect_all_present: bool,
    pub json: Option<String>,
    pub check: bool,
    pub reference: String,
    pub page: Option<String>,
    pub reproductions: Option<String>,
    pub list: bool,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            image: None,
            ensembl_vep: None,
            vep: None,
            vep_names_a_directory: false,
            perl_cache_dir: Vec::new(),
            plugins_dir: None,
            tool_dirs: Vec::new(),
            extra_dir: Vec::new(),
            entry: Vec::new(),
            jobs: 1,
            work_dir: None,
            keep: false,
            timeout: 1800,
            expect_all_present: false,
            json: None,
            check: false,
            reference: "116.3".to_string(),
            page: None,
            reproductions: None,
            list: false,
        }
    }
}

/// What a command line asks for.
#[derive(Debug, PartialEq)]
pub enum Command {
    Run(Box<Args>),
    Help,
}

/// The long options in the order a usage error lists them, each with whether it takes a
/// value.
const LONG_OPTIONS: [(&str, bool); 20] = [
    ("--help", false),
    ("--image", true),
    ("--ensembl-vep", true),
    ("--vep", true),
    ("--perl-cache-dir", true),
    ("--plugins-dir", true),
    ("--tool-dir", true),
    ("--extra-dir", true),
    ("--entry", true),
    ("--jobs", true),
    ("--work-dir", true),
    ("--keep", false),
    ("--timeout", true),
    ("--expect-all-present", false),
    ("--json", true),
    ("--check", false),
    ("--reference", true),
    ("--page", true),
    ("--reproductions", true),
    ("--list", false),
];

/// How one argument reads before any is consumed.
enum Reading<'a> {
    /// A value for the option before it, or a stray positional.
    Argument,
    /// One option, whether it takes a value, and the text after `=` when the argument
    /// carries one.
    Option(&'static str, bool, Option<&'a str>),
    /// `-h` with more characters behind it.
    ShortHelp(&'a str),
    /// A prefix that more than one option starts with.
    Ambiguous(Vec<&'static str>),
    /// Shaped like an option and naming none.
    Unknown,
}

/// Reads one argument other than the `--` marker. An argument that starts with `-` is an
/// option unless it is `-` alone, starts like a negative number (`-` or `-.` and a digit)
/// or holds a space; such a one is a value.
fn read(arg: &str) -> Reading<'_> {
    if !arg.starts_with('-') || arg.len() == 1 {
        return Reading::Argument;
    }
    let (head, inline) = match arg.split_once('=') {
        Some((head, inline)) => (head, Some(inline)),
        None => (arg, None),
    };
    if head == "-h" {
        return Reading::Option("--help", false, inline);
    }
    if arg.starts_with("--") {
        let exact = LONG_OPTIONS.iter().find(|(option, _)| *option == head);
        let matches: Vec<&(&'static str, bool)> = match exact {
            Some(option) => vec![option],
            None => LONG_OPTIONS
                .iter()
                .filter(|(option, _)| option.starts_with(head))
                .collect(),
        };
        match matches.as_slice() {
            [] => {}
            [(option, takes_value)] => return Reading::Option(option, *takes_value, inline),
            _ => return Reading::Ambiguous(matches.iter().map(|(option, _)| *option).collect()),
        }
    } else if let Some(rest) = arg.strip_prefix("-h") {
        return Reading::ShortHelp(rest);
    }
    let bytes = arg.as_bytes();
    let negative_number = matches!(
        bytes,
        [b'-', b'0'..=b'9', ..] | [b'-', b'.', b'0'..=b'9', ..]
    );
    if negative_number || arg.contains(' ') {
        Reading::Argument
    } else {
        Reading::Unknown
    }
}

/// The characters behind a `-h`. The option takes no value, so they read as more one-letter
/// options: a further `h` is help again, a letter that names no option is left unread and
/// the help is shown all the same, and a `-` or an `=` there is a value the option refuses.
fn short_help(mut rest: &str) -> Result<(), String> {
    let refused = |value: &str| format!("argument -h/--help: ignored explicit argument '{value}'");
    loop {
        if rest.starts_with('-') {
            return Err(refused(rest));
        }
        if let Some(value) = rest.strip_prefix("h=") {
            return Err(refused(value));
        }
        match rest.strip_prefix('h') {
            Some(more) if !more.is_empty() => rest = more,
            _ => return Ok(()),
        }
    }
}

fn int(name: &str, text: &str) -> Result<i64, String> {
    py_int_clamped(text).ok_or_else(|| format!("argument {name}: invalid int value: '{text}'"))
}

pub fn parse(argv: &[String]) -> Result<Command, String> {
    let mut args = Args::default();
    let mut unread: Vec<String> = Vec::new();
    let mut index = 0;
    while index < argv.len() {
        let arg = argv[index].as_str();
        index += 1;
        if arg == "--" {
            // The marker makes every later argument a positional, and the parser takes none.
            unread.extend(argv[index - 1..].iter().cloned());
            break;
        }
        let (name, takes_value, inline) = match read(arg) {
            Reading::Argument | Reading::Unknown => {
                unread.push(arg.to_string());
                continue;
            }
            Reading::Ambiguous(matches) => {
                return Err(format!(
                    "ambiguous option: {arg} could match {}",
                    matches.join(", ")
                ))
            }
            Reading::ShortHelp(rest) => return short_help(rest).map(|()| Command::Help),
            Reading::Option(name, takes_value, inline) => (name, takes_value, inline),
        };
        if !takes_value {
            if let Some(value) = inline {
                let shown = if name == "--help" { "-h/--help" } else { name };
                return Err(format!(
                    "argument {shown}: ignored explicit argument '{value}'"
                ));
            }
            match name {
                "--help" => return Ok(Command::Help),
                "--keep" => args.keep = true,
                "--expect-all-present" => args.expect_all_present = true,
                "--check" => args.check = true,
                "--list" => args.list = true,
                _ => unreachable!("a flag of LONG_OPTIONS"),
            }
            continue;
        }
        let value = match inline {
            Some(value) => value,
            None => match argv.get(index).map(String::as_str) {
                Some(next) if next != "--" && matches!(read(next), Reading::Argument) => {
                    index += 1;
                    next
                }
                _ => return Err(format!("argument {name}: expected one argument")),
            },
        };
        let path = || pypath::normalise(value);
        match name {
            "--image" => args.image = Some(value.to_string()),
            "--ensembl-vep" => args.ensembl_vep = Some(value.to_string()),
            "--vep" => {
                args.vep = Some(path());
                args.vep_names_a_directory = value.contains('/');
            }
            "--perl-cache-dir" => args.perl_cache_dir.push(path()),
            "--plugins-dir" => args.plugins_dir = Some(path()),
            "--tool-dir" => args.tool_dirs.push(path()),
            "--extra-dir" => args.extra_dir.push(value.to_string()),
            "--entry" => args.entry.push(value.to_string()),
            "--jobs" => args.jobs = int(name, value)?,
            "--work-dir" => args.work_dir = Some(path()),
            "--timeout" => args.timeout = int(name, value)?,
            "--json" => args.json = Some(path()),
            "--reference" => args.reference = value.to_string(),
            "--page" => args.page = Some(path()),
            "--reproductions" => args.reproductions = Some(path()),
            _ => unreachable!("an option of LONG_OPTIONS"),
        }
    }
    if !unread.is_empty() {
        return Err(format!("unrecognized arguments: {}", unread.join(" ")));
    }
    Ok(Command::Run(Box::new(args)))
}

#[cfg(test)]
#[path = "cli_tests.rs"]
mod tests;
