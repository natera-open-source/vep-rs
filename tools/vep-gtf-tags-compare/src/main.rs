// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-code agreement between the transcript tags of a GTF and the transcript attributes of
//! a JSON cache.
//!
//! `vep-cache-builder --gtf` turns the `tag` values on each GTF `transcript` line into the
//! four cache attributes VEP keeps (`basic` or `gencode_basic` to `gencode_basic`,
//! `gencode_primary`, `cds_start_NF`, `cds_end_NF`); a cache converted from Ensembl VEP's
//! holds the same four from the Ensembl core database. This binary reads a GTF the way the
//! builder does and the transcript shards of a JSON cache the way the runtime does, the GTF
//! and the shards side by side across the host's cores, and reports per code how many
//! transcripts agree on its presence and which disagree.

mod cache;
mod gtf;
mod input;
mod json;
mod pyfmt;
mod report;
mod threads;

use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

const USAGE: &str =
    "usage: vep-gtf-tags-compare [-h] --gtf GTF --cache CACHE [--out OUT] [--max-mismatches MAX_MISMATCHES]";
const HELP: &str = "\
Compare the transcript attributes a GTF yields with those a JSON cache carries.

Exit status is 0 when every code agrees on every transcript present on both sides,
1 otherwise, 2 on a usage error. Transcripts present on only one side are counted,
not treated as disagreements.

options:
  -h, --help            show this help message and exit
  --gtf GTF             GTF file, gzipped when its name ends in .gz
  --cache CACHE         JSON cache directory (holds transcripts/<chr>/)
  --out OUT             write the JSON report here as well as printing a summary
  --max-mismatches MAX_MISMATCHES
                        examples kept per code (default 50)

environment:
  VEP_GTF_TAGS_COMPARE_THREADS  worker threads (default: the host's available cores)";

#[derive(Debug, PartialEq)]
struct Args {
    gtf: PathBuf,
    cache: PathBuf,
    out: Option<PathBuf>,
    max_mismatches: i64,
}

#[derive(Debug, PartialEq)]
enum Parsed {
    Run(Args),
    Help,
    Usage(String),
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Opt {
    Help,
    Gtf,
    Cache,
    Out,
    MaxMismatches,
}

const LONG_OPTIONS: [(&str, Opt); 5] = [
    ("--help", Opt::Help),
    ("--gtf", Opt::Gtf),
    ("--cache", Opt::Cache),
    ("--out", Opt::Out),
    ("--max-mismatches", Opt::MaxMismatches),
];

/// What one argument is before it is used, by argparse's rules.
#[derive(Debug, PartialEq)]
enum Token<'a> {
    /// A value: an option takes it, and on its own it is an argument too many.
    Value,
    /// An option, with the text behind its `=` when it has one.
    Option(Opt, Option<&'a [u8]>),
    /// `-h` with letters behind it, which argparse reads as further one-letter options.
    HelpCluster(&'a [u8]),
    /// `--=...`, a prefix of every long option, or `-=...`, a prefix of every option.
    Ambiguous,
    /// A dash-led argument that names no option.
    Unknown,
}

/// Classifies one argument. An option is one of the five names, an unambiguous prefix of a
/// long one (the names share no first letter), or either with `=value` behind it. An
/// argument that is empty, `-`, dash-led like a negative number (`-5`, `-.5`) or holds a
/// space is a value.
fn classify(arg: &[u8]) -> Token<'_> {
    if arg.first() != Some(&b'-') || arg.len() == 1 {
        return Token::Value;
    }
    if arg == b"-h" {
        return Token::Option(Opt::Help, None);
    }
    let (name, explicit) = match memchr::memchr(b'=', arg) {
        Some(i) => (&arg[..i], Some(&arg[i + 1..])),
        None => (arg, None),
    };
    if name == b"-h" {
        return Token::Option(Opt::Help, explicit);
    }
    if name == b"-" {
        return Token::Ambiguous;
    }
    if arg.starts_with(b"--") {
        if name == b"--" {
            return Token::Ambiguous;
        }
        let named = LONG_OPTIONS
            .iter()
            .find(|(long, _)| long.as_bytes().starts_with(name));
        if let Some((_, opt)) = named {
            return Token::Option(*opt, explicit);
        }
    } else if arg.starts_with(b"-h") {
        return Token::HelpCluster(&arg[2..]);
    }
    let digit = |i: usize| arg.get(i).is_some_and(u8::is_ascii_digit);
    if digit(1) || (arg[1] == b'.' && digit(2)) || arg.contains(&b' ') {
        return Token::Value;
    }
    Token::Unknown
}

/// Reads the command line as argparse reads it, left to right: `--help` ends the reading
/// where it stands, a missing or unreadable value is an error where it stands, and an
/// argument too many or a missing required option is an error at the end.
fn parse_args(argv: &[OsString]) -> Parsed {
    let mut gtf = None;
    let mut cache = None;
    let mut out = None;
    let mut max_mismatches = 50;
    let mut surplus: Option<&OsStr> = None;
    let mut rest = argv.iter().map(OsString::as_os_str);
    while let Some(arg) = rest.next() {
        if arg == "--" {
            // Everything from here on is a value, and no option is left to take one.
            surplus = surplus.or(Some(arg));
            break;
        }
        let (opt, explicit) = match classify(arg.as_bytes()) {
            Token::Option(opt, explicit) => (opt, explicit),
            Token::HelpCluster(letters) if letters[0] != b'-' => return Parsed::Help,
            Token::HelpCluster(letters) => {
                return Parsed::Usage(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    String::from_utf8_lossy(letters)
                ));
            }
            Token::Ambiguous => {
                // `-` before the `=` is a prefix of `-h` as well.
                let short = if arg.as_bytes().starts_with(b"--") {
                    ""
                } else {
                    "-h, "
                };
                return Parsed::Usage(format!(
                    "ambiguous option: {} could match {short}--help, --gtf, --cache, --out, --max-mismatches",
                    arg.to_string_lossy()
                ));
            }
            Token::Value | Token::Unknown => {
                surplus = surplus.or(Some(arg));
                continue;
            }
        };
        let name = LONG_OPTIONS
            .iter()
            .find(|(_, o)| *o == opt)
            .map_or("", |(long, _)| long);
        if opt == Opt::Help {
            return match explicit {
                None => Parsed::Help,
                Some(text) => Parsed::Usage(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    String::from_utf8_lossy(text)
                )),
            };
        }
        let value = match explicit {
            Some(text) => OsStr::from_bytes(text),
            None => {
                let next = rest.clone().next();
                match next {
                    Some(next) if next != "--" && classify(next.as_bytes()) == Token::Value => {
                        rest.next();
                        next
                    }
                    _ => return Parsed::Usage(format!("argument {name}: expected one argument")),
                }
            }
        };
        match opt {
            Opt::Gtf => gtf = Some(PathBuf::from(value)),
            Opt::Cache => cache = Some(PathBuf::from(value)),
            Opt::Out => out = Some(PathBuf::from(value)),
            Opt::MaxMismatches => match value.to_str().and_then(pyfmt::py_int) {
                Some(bound) => max_mismatches = bound,
                None => {
                    return Parsed::Usage(format!(
                        "argument --max-mismatches: invalid int value: '{}'",
                        value.to_string_lossy()
                    ));
                }
            },
            Opt::Help => unreachable!("--help is handled above"),
        }
    }
    let (gtf, cache) = match (gtf, cache) {
        (Some(gtf), Some(cache)) => (gtf, cache),
        (gtf, cache) => {
            let missing: Vec<&str> = [(gtf, "--gtf"), (cache, "--cache")]
                .iter()
                .filter(|(value, _)| value.is_none())
                .map(|(_, name)| *name)
                .collect();
            return Parsed::Usage(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ));
        }
    };
    if let Some(arg) = surplus {
        return Parsed::Usage(format!("unrecognized arguments: {}", arg.to_string_lossy()));
    }
    Parsed::Run(Args {
        gtf,
        cache,
        // An empty `--out` writes no report.
        out: out.filter(|path| !path.as_os_str().is_empty()),
        max_mismatches,
    })
}

/// What one reading task returns: the GTF's code sets, or the transcripts of one shard.
enum Side {
    Gtf(gtf::CodeSets),
    Shard(cache::Transcripts),
}

fn read_inputs(args: &Args) -> Result<(gtf::CodeSets, gtf::CodeSets), String> {
    let threads = threads::thread_count()?;
    let shards = cache::list_shards(&args.cache);
    let named = |path: &Path, e: io::Error| format!("{}: {e}", path.display());
    // Task 0 reads the GTF and task `i` shard `i - 1`, so the one long stream and the many
    // shards share the threads.
    let sides = threads::run_indexed(threads, shards.len() + 1, |task| {
        match task.checked_sub(1) {
            None => input::open(&args.gtf)
                .and_then(gtf::read_gtf)
                .map(Side::Gtf)
                .map_err(|e| named(&args.gtf, e)),
            Some(shard) => cache::read_shard(&shards[shard])
                .map(Side::Shard)
                .map_err(|e| named(&shards[shard], e)),
        }
    })?;
    let mut from_gtf = gtf::CodeSets::new();
    let mut from_shards = Vec::with_capacity(shards.len());
    for side in sides {
        match side {
            Side::Gtf(sets) => from_gtf = sets,
            Side::Shard(transcripts) => from_shards.push(transcripts),
        }
    }
    Ok((from_gtf, cache::merge(&from_shards, &shards)?))
}

/// Runs the comparison: the stdout lines first, then the report. True on exact agreement.
fn run(args: &Args) -> Result<bool, String> {
    let (from_gtf, from_cache) = read_inputs(args)?;
    let report = report::compare(&from_gtf, &from_cache, args.max_mismatches);
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(&report::summary(&report))
        .and_then(|()| stdout.flush())
        .map_err(|e| format!("stdout: {e}"))?;
    if let Some(path) = &args.out {
        let named = |e: io::Error| format!("{}: {e}", path.display());
        let mut out = BufWriter::new(File::create(path).map_err(named)?);
        report::write_json(&report, &mut out)
            .and_then(|()| out.flush())
            .map_err(named)?;
    }
    let exact = report.exact_agreement();
    drop(report);
    // The process ends behind this run: freeing the ids one by one would only delay it.
    std::mem::forget((from_gtf, from_cache));
    Ok(exact)
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match parse_args(&argv) {
        Parsed::Run(args) => args,
        Parsed::Help => {
            println!("{USAGE}\n\n{HELP}");
            return ExitCode::SUCCESS;
        }
        Parsed::Usage(e) => {
            eprintln!("{USAGE}\nvep-gtf-tags-compare: error: {e}");
            return ExitCode::from(2);
        }
    };
    if !args.gtf.exists() || !args.cache.is_dir() {
        eprintln!(
            "ERROR: [compare_gtf_tags_to_cache] --gtf must be a file and --cache a directory"
        );
        return ExitCode::from(2);
    }
    match run(&args) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("ERROR: [vep-gtf-tags-compare] {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Parsed {
        parse_args(&args.iter().map(OsString::from).collect::<Vec<_>>())
    }

    fn run_with(gtf: &str, cache: &str, out: Option<&str>, max_mismatches: i64) -> Parsed {
        Parsed::Run(Args {
            gtf: PathBuf::from(gtf),
            cache: PathBuf::from(cache),
            out: out.map(PathBuf::from),
            max_mismatches,
        })
    }

    #[test]
    fn options_take_both_spellings_and_default_the_bound_to_50() {
        assert_eq!(
            parse(&["--gtf", "g", "--cache", "c"]),
            run_with("g", "c", None, 50)
        );
        assert_eq!(
            parse(&[
                "--cache=c",
                "--out",
                "o",
                "--gtf=g",
                "--max-mismatches",
                "7"
            ]),
            run_with("g", "c", Some("o"), 7)
        );
        assert_eq!(
            parse(&["--gtf", "first", "--gtf", "g", "--cache", "c", "--out="]),
            run_with("g", "c", None, 50)
        );
    }

    #[test]
    fn unambiguous_prefixes_name_their_option() {
        assert_eq!(
            parse(&["--g", "g", "--c=c", "--o", "o", "--m", "3"]),
            run_with("g", "c", Some("o"), 3)
        );
        assert_eq!(
            parse(&["--gt=g", "--cach", "c", "--max-mis=0"]),
            run_with("g", "c", None, 0)
        );
        assert_eq!(parse(&["--h"]), Parsed::Help);
    }

    #[test]
    fn values_follow_argparse_s_reading_of_dash_led_arguments() {
        assert_eq!(
            parse(&["--gtf", "-", "--cache", "-5", "--max-mismatches", "-1"]),
            run_with("-", "-5", None, -1)
        );
        assert_eq!(
            parse(&["--gtf", "-my file", "--cache", "-.5x", "--out", ""]),
            run_with("-my file", "-.5x", None, 50)
        );
        assert_eq!(
            parse(&["--gtf=--cache", "--cache=-x"]),
            run_with("--cache", "-x", None, 50)
        );
        for argv in [
            &["--gtf", "-x", "--cache", "c"][..],
            &["--gtf", "--cache", "c"],
            &["--gtf", "g", "--cache"],
            &["--gtf", "--", "g", "--cache", "c"],
            &["--gtf", "-.x", "--cache", "c"],
        ] {
            assert!(matches!(parse(argv), Parsed::Usage(_)), "{argv:?}");
        }
    }

    #[test]
    fn the_bound_is_read_as_pythons_int() {
        let bound = |value: &str| parse(&["--gtf", "g", "--cache", "c", "--max-mismatches", value]);
        assert_eq!(bound(" 1_0 "), run_with("g", "c", None, 10));
        assert_eq!(bound("+0"), run_with("g", "c", None, 0));
        assert_eq!(
            bound("99999999999999999999999999"),
            run_with("g", "c", None, i64::MAX)
        );
        for value in ["", "x", "5.0", "1e3", "0x10", "-.5", "-5x"] {
            assert!(matches!(bound(value), Parsed::Usage(_)), "{value:?}");
        }
    }

    #[test]
    fn help_ends_the_reading_where_it_stands() {
        assert_eq!(parse(&["-h"]), Parsed::Help);
        assert_eq!(parse(&["--help"]), Parsed::Help);
        assert_eq!(parse(&["--bogus", "extra", "--help"]), Parsed::Help);
        assert_eq!(
            parse(&["--gtf", "g", "-h", "--max-mismatches", "x"]),
            Parsed::Help
        );
        assert_eq!(parse(&["-hx"]), Parsed::Help);
        assert_eq!(parse(&["-hh"]), Parsed::Help);
        for argv in [
            &["--max-mismatches", "x", "--help"][..],
            &["--gtf", "--help"],
            &["--", "--help"],
            &["--help=x"],
            &["--help="],
            &["-h=x"],
            &["-h="],
            &["-h-x"],
            &["--=x", "--help"],
        ] {
            assert!(matches!(parse(argv), Parsed::Usage(_)), "{argv:?}");
        }
    }

    #[test]
    fn a_dash_before_an_equals_sign_is_a_prefix_of_every_option() {
        let ambiguous = |arg: &str, short: &str| {
            Parsed::Usage(format!(
                "ambiguous option: {arg} could match {short}--help, --gtf, --cache, --out, --max-mismatches"
            ))
        };
        assert_eq!(parse(&["-=x", "--help"]), ambiguous("-=x", "-h, "));
        assert_eq!(
            parse(&["--gtf", "g", "--cache", "c", "-="]),
            ambiguous("-=", "-h, ")
        );
        assert_eq!(parse(&["--=x", "--help"]), ambiguous("--=x", ""));
        // What is read before it is reported first, and it is no value for an option.
        assert_eq!(parse(&["-h", "-=x"]), Parsed::Help);
        assert_eq!(
            parse(&["--gtf", "-=x"]),
            Parsed::Usage("argument --gtf: expected one argument".to_string())
        );
    }

    #[test]
    fn missing_options_and_arguments_too_many_are_usage_errors() {
        for argv in [
            &[][..],
            &["--gtf", "g"],
            &["--cache", "c"],
            &["--gtf", "g", "--cache", "c", "extra"],
            &["extra", "--gtf", "g", "--cache", "c"],
            &["--gtf", "g", "--cache", "c", "--bogus"],
            &["--gtf", "g", "--cache", "c", "--bogus=1"],
            &["--gtf", "g", "--cache", "c", "-x"],
            &["--gtf", "g", "--cache", "c", "-gtf", "g"],
            &["--gtf", "g", "--cache", "c", "--"],
            &["--gtf", "g", "--", "--cache", "c"],
            &["--gtf", "g", "--cache", "c", "--outx", "o"],
            &["--gtf", "g", "--cache", "c", ""],
            &["--gtf", "g", "--cache", "c", "-"],
        ] {
            assert!(matches!(parse(argv), Parsed::Usage(_)), "{argv:?}");
        }
    }

    #[test]
    fn values_that_are_not_utf8_reach_their_option_as_given() {
        let raw = OsStr::from_bytes(b"dir\xff/genes.gtf");
        let mut inline = OsString::from("--gtf=");
        inline.push(raw);
        for argv in [
            vec![inline, "--cache".into(), "c".into()],
            vec![
                "--gtf".into(),
                raw.to_os_string(),
                "--cache".into(),
                "c".into(),
            ],
        ] {
            match parse_args(&argv) {
                Parsed::Run(args) => assert_eq!(args.gtf.as_os_str(), raw),
                other => panic!("{other:?}"),
            }
        }
    }
}
