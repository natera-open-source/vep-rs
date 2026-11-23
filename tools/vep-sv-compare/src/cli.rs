// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line: the options with their defaults, the help, the usage errors and the
//! exit codes.
//!
//! Exit codes: 2 on a usage error, 1 when the comparison refuses its inputs or fails, 0
//! otherwise, the help included; 120 in place of any of them when stdout or stderr cannot
//! take what the run wrote ([`streams`]).

pub mod streams;

use std::ffi::{OsStr, OsString};
use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::LazyLock;

use regex::Regex;

use self::streams::Streams;
use crate::compare::{run, Stop};
use crate::masks::{DEFAULT_REFERENCE_RELEASE, MAX_SV_SIZE, REFERENCE_RELEASES};
use crate::text::{py_path, py_repr};
use crate::threads::thread_count;

/// The program's name in its own messages.
pub const PROG: &str = "vep-sv-compare";

/// The help: the text `argparse` writes for these options under the program's name, 80
/// columns wide and without colour.
pub const HELP: &str = include_str!("help.txt");

/// The usage block that opens the help and every usage error: the help's first paragraph.
pub fn usage() -> &'static str {
    let end = HELP
        .find("\n\n")
        .expect("a blank line ends the usage block");
    &HELP[..=end]
}

/// The assemblies `--assembly` takes, the first the default.
pub const ASSEMBLIES: [&str; 2] = ["grch37", "grch38"];

/// The engine `--scored-engine` names when it is not given.
pub const DEFAULT_SCORED_ENGINE: &str = "vep-rs";

/// What a comparison is run with.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub assembly: String,
    /// The directory of the input VCF files.
    pub input_dir: PathBuf,
    /// The directory of the reference's outputs.
    pub perl_dir: PathBuf,
    /// The directory of the scored engine's outputs.
    pub rust_dir: PathBuf,
    /// The directory the reports are written to.
    pub output_dir: PathBuf,
    /// The vep-rs JSON cache directory as it was written on the command line; `None` when
    /// the option is absent or empty.
    pub vep_rs_cache: Option<OsString>,
    pub scored_engine: String,
    pub reference_release: String,
    /// The `--max_sv_size` the reference ran with; -1 is no cap.
    pub reference_max_sv_size: i128,
    /// The worker-thread bound, from the environment and not the command line.
    pub threads: usize,
}

impl Options {
    /// The options of a comparison of `assembly` with every other option at its default.
    pub fn for_assembly(assembly: &str) -> Options {
        let directory = |root: &str| PathBuf::from(format!("{root}/{assembly}"));
        Options {
            assembly: assembly.to_string(),
            input_dir: directory("tests/sv_validation"),
            perl_dir: directory("tmp/sv_perl"),
            rust_dir: directory("tmp/sv_rust"),
            output_dir: directory("tmp/sv_concordance"),
            vep_rs_cache: None,
            scored_engine: DEFAULT_SCORED_ENGINE.to_string(),
            reference_release: DEFAULT_REFERENCE_RELEASE.to_string(),
            reference_max_sv_size: MAX_SV_SIZE,
            threads: 1,
        }
    }
}

/// The option strings in the order a message lists them: the two of the help, then the
/// nine options that take a value.
const OPTION_STRINGS: [&str; 11] = [
    "-h",
    "--help",
    "--assembly",
    "--input-dir",
    "--perl-dir",
    "--rust-dir",
    "--output-dir",
    "--vep-rs-cache",
    "--scored-engine",
    "--reference-release",
    "--reference-max-sv-size",
];

/// The help option as an error names it.
const HELP_OPTION: &str = "-h/--help";

/// What a command line asks for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Request {
    /// A comparison run with these options.
    Compare(Box<Options>),
    /// The help.
    Help,
}

/// One reading of an argument as an option.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Reading<'a> {
    /// The option string the argument spells or is a prefix of.
    option: &'static str,
    /// Whether an `=` stands between the option's name and the text after it.
    equals: bool,
    /// The text after the option's name in the same argument.
    explicit: Option<&'a [u8]>,
}

/// What an argument is before any option is taken.
#[derive(Debug, PartialEq, Eq)]
enum Argument<'a> {
    /// No option: the value of an option that stands before it, else unrecognised.
    Positional,
    /// An argument that begins with `-` and names no option: unrecognised, and never a
    /// value.
    Unknown,
    /// The options the argument can name; more than one is an ambiguous prefix.
    Named(Vec<Reading<'a>>),
}

/// Reads one argument as `argparse` does: an option string whole, with `=value`, or as
/// a prefix of the text before `=`; `-h` also with further text directly after it. An
/// argument that names no option is positional when it does not begin with `-`, is `-`
/// alone, begins as a negative number or holds a blank.
fn read_argument(argument: &[u8]) -> Argument<'_> {
    static NEGATIVE_NUMBER: LazyLock<regex::bytes::Regex> =
        LazyLock::new(|| regex::bytes::Regex::new(r"^-\.?\d").expect("regex"));
    if !argument.starts_with(b"-") || argument.len() == 1 {
        return Argument::Positional;
    }
    let spelled = |name: &[u8]| {
        OPTION_STRINGS
            .iter()
            .copied()
            .find(|o| o.as_bytes() == name)
    };
    if let Some(option) = spelled(argument) {
        return Argument::Named(vec![Reading {
            option,
            equals: false,
            explicit: None,
        }]);
    }
    let (name, explicit) = match memchr::memchr(b'=', argument) {
        Some(i) => (&argument[..i], Some(&argument[i + 1..])),
        None => (argument, None),
    };
    if let (Some(option), Some(_)) = (spelled(name), explicit) {
        return Argument::Named(vec![Reading {
            option,
            equals: true,
            explicit,
        }]);
    }
    let one_dash = argument[1] != b'-';
    let readings: Vec<Reading> = OPTION_STRINGS
        .iter()
        .copied()
        .filter_map(|option| {
            if one_dash && option.as_bytes() == &argument[..2] {
                Some(Reading {
                    option,
                    equals: false,
                    explicit: Some(&argument[2..]),
                })
            } else if option.as_bytes().starts_with(name) {
                Some(Reading {
                    option,
                    equals: explicit.is_some(),
                    explicit,
                })
            } else {
                None
            }
        })
        .collect();
    if !readings.is_empty() {
        return Argument::Named(readings);
    }
    if NEGATIVE_NUMBER.is_match(argument) || argument.contains(&b' ') {
        return Argument::Positional;
    }
    Argument::Unknown
}

/// Whether an argument that names the help option asks for the help; the message of the
/// usage error it is otherwise. `-h` runs on into further one-letter options of the same
/// argument: another `h` is the help again, and any other letter is an unrecognised
/// option the help is printed before.
fn help_asked(reading: &Reading) -> Result<(), String> {
    let ignored = |explicit: &[u8]| {
        format!(
            "argument {HELP_OPTION}: ignored explicit argument {}",
            py_repr(explicit)
        )
    };
    let Some(mut explicit) = reading.explicit else {
        return Ok(());
    };
    if reading.option != "-h" {
        return Err(ignored(explicit));
    }
    let mut equals = reading.equals;
    loop {
        if equals || explicit.is_empty() || explicit[0] == b'-' {
            return Err(ignored(explicit));
        }
        if explicit[0] != b'h' || explicit.len() == 1 {
            return Ok(());
        }
        explicit = &explicit[1..];
        if explicit[0] == b'=' {
            equals = true;
            explicit = &explicit[1..];
        }
    }
}

/// An argument as Python writes it to stderr: a byte that is not UTF-8 is the escape of
/// the lone surrogate Python holds it as (`\udcNN`).
fn py_stderr_text(argument: &[u8]) -> String {
    let mut out = String::with_capacity(argument.len());
    for chunk in argument.utf8_chunks() {
        out.push_str(chunk.valid());
        for byte in chunk.invalid() {
            out.push_str(&format!("\\udc{byte:02x}"));
        }
    }
    out
}

/// The value of a decimal digit of any script.
fn decimal_digit(c: char) -> Option<u32> {
    static DIGIT: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\d$").expect("regex"));
    let is_digit = |c: char| {
        let mut buffer = [0u8; 4];
        DIGIT.is_match(c.encode_utf8(&mut buffer))
    };
    if let Some(digit) = c.to_digit(10) {
        return Some(digit);
    }
    if !is_digit(c) {
        return None;
    }
    // The decimal digits of a script are ten consecutive characters from its zero, and the
    // runs of neighbouring scripts follow each other in whole tens.
    let mut zero = c as u32;
    while zero > 0 && char::from_u32(zero - 1).is_some_and(is_digit) {
        zero -= 1;
    }
    Some((c as u32 - zero) % 10)
}

/// Python's `int()` of an argument: blanks around the number dropped, an optional sign,
/// decimal digits of any script with single underscores between them. `None` where Python
/// raises, and for a value outside 128 bits, which Python would take.
fn py_int(text: &str) -> Option<i128> {
    let ascii: String = text
        .chars()
        .map(|c| match decimal_digit(c) {
            Some(digit) if !c.is_ascii() => char::from_digit(digit, 10).expect("a digit"),
            _ => c,
        })
        .collect();
    crate::rules::py_int(&ascii)
}

/// What the arguments after the program's name ask for, or the message of the usage
/// error they are. The arguments are taken in their order: the help is asked for where
/// its option is reached, whatever follows it, and a defective option before it is the
/// usage error. An argument that names no option is reported once every option is read.
pub fn parse_args(arguments: Vec<OsString>) -> Result<Request, String> {
    let mut assembly = ASSEMBLIES[0].to_string();
    let mut directories: [Option<PathBuf>; 4] = [None, None, None, None];
    let mut vep_rs_cache: Option<OsString> = None;
    let mut scored_engine = DEFAULT_SCORED_ENGINE.to_string();
    let mut reference_release = DEFAULT_REFERENCE_RELEASE.to_string();
    let mut reference_max_sv_size = MAX_SV_SIZE;
    let mut unrecognized: Vec<&OsStr> = Vec::new();
    // Every argument from the first `--` on is positional, and none is taken.
    let options_end = arguments
        .iter()
        .position(|a| a.as_bytes() == b"--")
        .unwrap_or(arguments.len());
    let mut index = 0;
    while index < options_end {
        let argument = &arguments[index];
        index += 1;
        let readings = match read_argument(argument.as_bytes()) {
            Argument::Named(readings) => readings,
            Argument::Positional | Argument::Unknown => {
                unrecognized.push(argument);
                continue;
            }
        };
        let [reading] = readings.as_slice() else {
            let matches: Vec<&str> = readings.iter().map(|r| r.option).collect();
            return Err(format!(
                "ambiguous option: {} could match {}",
                py_stderr_text(argument.as_bytes()),
                matches.join(", ")
            ));
        };
        let option = reading.option;
        if option == "-h" || option == "--help" {
            return help_asked(reading).map(|()| Request::Help);
        }
        let value: &OsStr = match reading.explicit {
            Some(explicit) => OsStr::from_bytes(explicit),
            None => {
                let next = arguments[..options_end]
                    .get(index)
                    .filter(|next| read_argument(next.as_bytes()) == Argument::Positional);
                index += 1;
                next.ok_or_else(|| format!("argument {option}: expected one argument"))?
            }
        };
        let text = value.to_string_lossy().into_owned();
        let choice = |choices: &[&str]| -> Result<String, String> {
            if choices.contains(&text.as_str()) {
                return Ok(text.clone());
            }
            let listed: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
            Err(format!(
                "argument {option}: invalid choice: {} (choose from {})",
                py_repr(value.as_bytes()),
                listed.join(", ")
            ))
        };
        match option {
            "--assembly" => assembly = choice(&ASSEMBLIES)?,
            "--input-dir" => directories[0] = Some(py_path(value)),
            "--perl-dir" => directories[1] = Some(py_path(value)),
            "--rust-dir" => directories[2] = Some(py_path(value)),
            "--output-dir" => directories[3] = Some(py_path(value)),
            "--vep-rs-cache" => vep_rs_cache = Some(value.to_os_string()),
            "--scored-engine" => scored_engine = text,
            "--reference-release" => reference_release = choice(&REFERENCE_RELEASES)?,
            "--reference-max-sv-size" => {
                reference_max_sv_size = py_int(&text).ok_or_else(|| {
                    format!(
                        "argument {option}: invalid int value: {}",
                        py_repr(value.as_bytes())
                    )
                })?;
            }
            _ => unreachable!("an option that takes a value"),
        }
    }
    unrecognized.extend(arguments[options_end..].iter().map(OsString::as_os_str));
    if !unrecognized.is_empty() {
        let listed: Vec<String> = unrecognized
            .iter()
            .map(|argument| py_stderr_text(argument.as_bytes()))
            .collect();
        return Err(format!("unrecognized arguments: {}", listed.join(" ")));
    }
    let defaults = Options::for_assembly(&assembly);
    let [input_dir, perl_dir, rust_dir, output_dir] = directories;
    Ok(Request::Compare(Box::new(Options {
        input_dir: input_dir.unwrap_or(defaults.input_dir),
        perl_dir: perl_dir.unwrap_or(defaults.perl_dir),
        rust_dir: rust_dir.unwrap_or(defaults.rust_dir),
        output_dir: output_dir.unwrap_or(defaults.output_dir),
        vep_rs_cache: vep_rs_cache.filter(|cache| !cache.is_empty()),
        scored_engine,
        reference_release,
        reference_max_sv_size,
        ..defaults
    })))
}

/// Runs what the process arguments ask for, the comparison or the help on stdout, and
/// ends as [`Streams::close`] says.
pub fn main() -> ExitCode {
    ignore_file_size_signal();
    let mut streams = Streams::of_process();
    let code = command(std::env::args_os().skip(1).collect(), &mut streams);
    ExitCode::from(streams.close(code))
}

/// Sets the signal of a write past the file-size limit to ignored, as the interpreter sets
/// it when it starts. The write then fails with the system's error (`File too large`), as a
/// write to a pipe whose reader has gone fails, and the run ends as it does on any write a
/// file does not take.
fn ignore_file_size_signal() {
    // SAFETY: the call is made once, at the start of `main`, before any thread of the run
    // starts. It sets the disposition of one signal to "ignored" and installs no handler,
    // so no code of this process ever runs in a signal's context.
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGXFSZ, libc::SIG_IGN);
    }
}

/// Runs what `arguments` ask for on `streams` and gives the exit code of a run whose
/// streams take everything: 0, 1 or 2 as the module says.
///
/// No write to a stream ends the run here. The help and a usage error are written whether
/// the stream takes them or not; the lines of a refusal stop at the first that stderr does
/// not take. A refusal writes its lines while stdout still holds the text of the run, and
/// a failure writes its line once that text is flushed.
pub fn command<O: Write, E: Write>(arguments: Vec<OsString>, streams: &mut Streams<O, E>) -> u8 {
    let mut options = match parse_args(arguments) {
        Ok(Request::Compare(options)) => *options,
        Ok(Request::Help) => {
            let _ = streams.out.write_text(HELP.as_bytes());
            return 0;
        }
        Err(message) => {
            let _ = streams.err.write_text(usage().as_bytes());
            let _ = streams
                .err
                .write_text(format!("{PROG}: error: {message}\n").as_bytes());
            return 2;
        }
    };
    options.threads = match thread_count() {
        Ok(threads) => threads,
        Err(failure) => {
            let _ = streams.err.print(&format!("ERROR: [{PROG}] {failure}"));
            return 1;
        }
    };
    match run(&options, &mut streams.out) {
        Ok(record) => {
            // The process ends here: the sets are left to it and not freed one by one.
            std::mem::forget(record);
            0
        }
        Err(Stop::Refused { lines, record }) => {
            std::mem::forget(record);
            for line in lines {
                if streams.err.print(&line).is_err() {
                    break;
                }
            }
            1
        }
        Err(Stop::Failed(failure)) => {
            let _ = streams.out.flush();
            let _ = streams.err.print(&format!("ERROR: [{PROG}] {failure}"));
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::streams::tests::{Sink, EPIPE};
    use super::streams::{Held, FLUSH_FAILED, MIN_BLOCK_BYTES};
    use super::*;

    /// Streams over two files that record what they take: a stdout that is no terminal,
    /// and stderr.
    fn recording() -> (Streams<Sink, Sink>, Sink, Sink) {
        let (out, err) = (Sink::default(), Sink::default());
        let streams = Streams {
            out: Held::new(out.clone(), MIN_BLOCK_BYTES, false),
            err: Held::new(err.clone(), MIN_BLOCK_BYTES, true),
        };
        (streams, out, err)
    }

    #[test]
    fn the_help_goes_to_stdout_and_a_usage_error_to_stderr() {
        let (mut streams, out, err) = recording();
        assert_eq!(command(vec!["--help".into()], &mut streams), 0);
        assert_eq!(out.text(), "");
        assert_eq!(streams.close(0), 0);
        assert_eq!((out.text().as_str(), err.text().as_str()), (HELP, ""));

        let (mut streams, out, err) = recording();
        let arguments = vec!["--assembly".into(), "bad".into()];
        assert_eq!(command(arguments, &mut streams), 2);
        assert_eq!(
            err.text(),
            format!(
                "{}{PROG}: error: argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')\n",
                usage()
            )
        );
        assert_eq!(streams.close(2), 2);
        assert_eq!(out.text(), "");
    }

    #[test]
    fn a_stream_that_takes_nothing_ends_neither_the_help_nor_a_usage_error() {
        let (mut streams, out, err) = recording();
        out.refuses.set(Some(EPIPE));
        assert_eq!(command(vec!["--help".into()], &mut streams), 0);
        assert_eq!(streams.close(0), FLUSH_FAILED);
        assert_eq!(
            err.text(),
            "Exception ignored while flushing sys.stdout:\nBrokenPipeError: [Errno 32] Broken pipe\n"
        );

        let (mut streams, _, err) = recording();
        err.refuses.set(Some(EPIPE));
        let arguments = vec!["--assembly".into(), "bad".into()];
        assert_eq!(command(arguments, &mut streams), 2);
        assert_eq!(streams.close(2), FLUSH_FAILED);
        assert!(err.writes.borrow().is_empty());
    }

    fn request(arguments: &[&str]) -> Result<Request, String> {
        parse_args(arguments.iter().map(OsString::from).collect())
    }

    fn parse(arguments: &[&str]) -> Result<Options, String> {
        request(arguments).map(|request| match request {
            Request::Compare(options) => *options,
            Request::Help => panic!("the help is asked for: {arguments:?}"),
        })
    }

    fn bytes(arguments: &[&[u8]]) -> Result<Request, String> {
        parse_args(
            arguments
                .iter()
                .map(|a| OsStr::from_bytes(a).to_os_string())
                .collect(),
        )
    }

    #[test]
    fn every_option_has_its_default() {
        let options = parse(&[]).unwrap();
        assert_eq!(options, Options::for_assembly("grch37"));
        assert_eq!(
            options.input_dir,
            PathBuf::from("tests/sv_validation/grch37")
        );
        assert_eq!(options.perl_dir, PathBuf::from("tmp/sv_perl/grch37"));
        assert_eq!(options.rust_dir, PathBuf::from("tmp/sv_rust/grch37"));
        assert_eq!(
            options.output_dir,
            PathBuf::from("tmp/sv_concordance/grch37")
        );
        assert_eq!(options.vep_rs_cache, None);
        assert_eq!(options.scored_engine, "vep-rs");
        assert_eq!(options.reference_release, "116.3");
        assert_eq!(options.reference_max_sv_size, 10_000_000);
        let options = parse(&["--assembly", "grch38", "--perl-dir", "p"]).unwrap();
        assert_eq!(
            options.input_dir,
            PathBuf::from("tests/sv_validation/grch38")
        );
        assert_eq!(options.perl_dir, PathBuf::from("p"));
        assert_eq!(
            options.output_dir,
            PathBuf::from("tmp/sv_concordance/grch38")
        );
    }

    #[test]
    fn an_option_is_spelled_whole_with_an_equals_sign_or_as_a_prefix() {
        let options = parse(&[
            "--input-dir=/data/in//x/",
            "--perl",
            "./perl",
            "--ru=rust",
            "--o",
            "",
            "--vep-rs-cache",
            "cache//dir/",
            "--scored-engine=fastVEP",
            "--reference-rel",
            "115.2",
            "--reference-m=-1",
            "--assembly=grch38",
        ])
        .unwrap();
        assert_eq!(options.assembly, "grch38");
        assert_eq!(options.input_dir, PathBuf::from("/data/in/x"));
        assert_eq!(options.perl_dir, PathBuf::from("perl"));
        assert_eq!(options.rust_dir, PathBuf::from("rust"));
        assert_eq!(options.output_dir, PathBuf::from("."));
        assert_eq!(options.vep_rs_cache, Some(OsString::from("cache//dir/")));
        assert_eq!(options.scored_engine, "fastVEP");
        assert_eq!(options.reference_release, "115.2");
        assert_eq!(options.reference_max_sv_size, -1);
        assert_eq!(parse(&["--vep-rs-cache", ""]).unwrap().vep_rs_cache, None);
        assert_eq!(parse(&["--vep-rs-cache="]).unwrap().vep_rs_cache, None);
    }

    #[test]
    fn the_last_of_a_repeated_option_wins() {
        let options = parse(&[
            "--perl-dir",
            "first",
            "--rust-dir",
            "x",
            "--perl-dir",
            "second",
            "--rust-di",
            "y",
            "--scored",
            "a",
            "--scored-engine",
            "b",
        ])
        .unwrap();
        assert_eq!(options.perl_dir, PathBuf::from("second"));
        assert_eq!(options.rust_dir, PathBuf::from("y"));
        assert_eq!(options.scored_engine, "b");
    }

    #[test]
    fn the_size_cap_is_an_integer_as_python_reads_one() {
        let cap =
            |text: &str| parse(&["--reference-max-sv-size", text]).map(|o| o.reference_max_sv_size);
        assert_eq!(cap("5"), Ok(5));
        assert_eq!(cap("-1"), Ok(-1));
        assert_eq!(cap(" 1_0 "), Ok(10));
        assert_eq!(cap("+7"), Ok(7));
        assert_eq!(cap("\u{ff11}\u{ff12}"), Ok(12));
        assert_eq!(cap("\u{663}\u{660}"), Ok(30));
        for refused in ["abc", "1.5", "", "1__0", "0x10", "\u{2167}"] {
            assert_eq!(
                cap(refused),
                Err(format!(
                    "argument --reference-max-sv-size: invalid int value: '{refused}'"
                )),
                "{refused}"
            );
        }
        assert_eq!(
            parse(&["--reference-max-sv-size=-1"])
                .unwrap()
                .reference_max_sv_size,
            -1
        );
    }

    #[test]
    fn usage_errors_carry_the_message_of_the_first_defect() {
        let message = |arguments: &[&str]| parse(arguments).unwrap_err();
        assert_eq!(
            message(&["--assembly", "x"]),
            "argument --assembly: invalid choice: 'x' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            message(&["--reference-release", "117"]),
            "argument --reference-release: invalid choice: '117' (choose from '115.2', '116.2', '116.3')"
        );
        assert_eq!(
            message(&["--assembly", "grch37", "--assembly", "bad"]),
            "argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            message(&["--re", "1"]),
            "ambiguous option: --re could match --reference-release, --reference-max-sv-size"
        );
        assert_eq!(
            message(&["--assembly", "grch38", "--r=1"]),
            "ambiguous option: --r=1 could match --rust-dir, --reference-release, --reference-max-sv-size"
        );
        assert_eq!(
            message(&["--perl-dir"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(
            message(&["--perl-dir", "--rust-dir", "x"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(
            message(&["--perl-dir", "-x"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(message(&["x"]), "unrecognized arguments: x");
        assert_eq!(
            message(&["x", "y", "--bogus"]),
            "unrecognized arguments: x y --bogus"
        );
        assert_eq!(message(&["-a"]), "unrecognized arguments: -a");
        assert_eq!(
            message(&["--input-dir", "/nonexistent", "extra"]),
            "unrecognized arguments: extra"
        );
        assert_eq!(
            message(&["extra", "--assembly", "bad"]),
            "argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            message(&["--perl-dir", "--assembly=x y"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(
            message(&["--perl-dir", "--", "x"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(
            message(&["--input-dir", "d", "--"]),
            "unrecognized arguments: --"
        );
    }

    /// The outcome of each argument list: `None` for the help, else the usage error.
    fn assert_outcomes(cases: &[(&[&str], Option<&str>)]) {
        for (arguments, expected) in cases {
            let outcome = match request(arguments) {
                Ok(Request::Help) => None,
                Ok(Request::Compare(_)) => panic!("neither help nor error: {arguments:?}"),
                Err(message) => Some(message),
            };
            assert_eq!(outcome.as_deref(), *expected, "{arguments:?}");
        }
    }

    #[test]
    fn the_help_is_asked_for_where_its_option_is_reached() {
        assert_outcomes(&[
            (&["--help"], None),
            (&["-h"], None),
            (&["--h"], None),
            (&["--he"], None),
            (&["--hel"], None),
            (&["x", "--help"], None),
            (&["-h", "x"], None),
            (&["-a", "-h"], None),
            (&["--unknown=1", "-h"], None),
            (&["", "-h"], None),
            (&["-", "-h"], None),
            (&["-1", "-h"], None),
            (&["--assembly", "grch38", "-h"], None),
            (&["--perl-dir", "", "-h"], None),
            (&["--i", "x", "--p", "y", "--h"], None),
            (&["-h", "--"], None),
            (&["-h", "--assembly", "bad"], None),
            (&["--help", "--re", "1"], None),
            (&["-h", "--r"], None),
            (&["-h", "--=x"], None),
            (&["-h", "-=x"], None),
            (&["-h", "--perl-dir"], None),
            (&["-h", "--reference-max-sv-size", "abc"], None),
        ]);
    }

    #[test]
    fn a_defect_before_the_help_option_is_the_usage_error() {
        assert_outcomes(&[
            (&["--assembly", "bad", "--help"], Some("argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')")),
            (&["--re", "1", "--help"], Some("ambiguous option: --re could match --reference-release, --reference-max-sv-size")),
            (&["--=x", "-h"], Some("ambiguous option: --=x could match --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size")),
            (&["-=x", "-h"], Some("ambiguous option: -=x could match -h, --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size")),
            (&["--reference-max-sv-size", "abc", "-h"], Some("argument --reference-max-sv-size: invalid int value: 'abc'")),
        ]);
    }

    #[test]
    fn the_help_option_is_no_value_and_is_positional_after_a_double_dash() {
        assert_outcomes(&[
            (
                &["--perl-dir", "-h"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "--help"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "--h"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "-hx"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "-h x"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "--he=x y"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "-=x"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--perl-dir", "-= x"],
                Some("argument --perl-dir: expected one argument"),
            ),
            (
                &["--reference-max-sv-size", "-h"],
                Some("argument --reference-max-sv-size: expected one argument"),
            ),
            (&["--", "--help"], Some("unrecognized arguments: -- --help")),
            (&["--", "x", "-h"], Some("unrecognized arguments: -- x -h")),
        ]);
        assert_eq!(
            parse(&["--perl-dir=-h"]).unwrap().perl_dir,
            PathBuf::from("-h")
        );
        assert_eq!(
            parse(&["--scored-engine=--help"]).unwrap().scored_engine,
            "--help"
        );
        assert_eq!(
            parse(&["--perl-dir", "--help x"]).unwrap().perl_dir,
            PathBuf::from("--help x")
        );
    }

    #[test]
    fn text_after_the_short_help_option_is_further_short_options() {
        assert_outcomes(&[
            (&["-hx"], None),
            (&["-hh"], None),
            (&["-ha"], None),
            (&["-h1"], None),
            (&["-help"], None),
            (&["-h x"], None),
            (&["-hx y"], None),
            (&["-hhhx"], None),
            (&["-hxh"], None),
            (&["-hx=1"], None),
            (&["-hu=1"], None),
            (&["-h\u{e9}"], None),
        ]);
    }

    #[test]
    fn the_help_option_takes_no_value() {
        assert_outcomes(&[
            (
                &["--help=x"],
                Some("argument -h/--help: ignored explicit argument 'x'"),
            ),
            (
                &["--h=x"],
                Some("argument -h/--help: ignored explicit argument 'x'"),
            ),
            (
                &["--help="],
                Some("argument -h/--help: ignored explicit argument ''"),
            ),
            (
                &["--he="],
                Some("argument -h/--help: ignored explicit argument ''"),
            ),
            (
                &["--he=x y"],
                Some("argument -h/--help: ignored explicit argument 'x y'"),
            ),
            (
                &["-h=x"],
                Some("argument -h/--help: ignored explicit argument 'x'"),
            ),
            (
                &["-h="],
                Some("argument -h/--help: ignored explicit argument ''"),
            ),
            (
                &["-h=x y"],
                Some("argument -h/--help: ignored explicit argument 'x y'"),
            ),
            (
                &["-h-"],
                Some("argument -h/--help: ignored explicit argument '-'"),
            ),
            (
                &["-h--"],
                Some("argument -h/--help: ignored explicit argument '--'"),
            ),
            (
                &["-h-h"],
                Some("argument -h/--help: ignored explicit argument '-h'"),
            ),
            (
                &["-hh=1"],
                Some("argument -h/--help: ignored explicit argument '1'"),
            ),
            (
                &["-hh="],
                Some("argument -h/--help: ignored explicit argument ''"),
            ),
        ]);
    }

    #[test]
    fn an_argument_that_only_resembles_the_help_option_is_unrecognised() {
        assert_outcomes(&[
            (&["--help x"], Some("unrecognized arguments: --help x")),
            (&["-ah"], Some("unrecognized arguments: -ah")),
            (&["-qh"], Some("unrecognized arguments: -qh")),
            (&["-H"], Some("unrecognized arguments: -H")),
            (&["--HELP"], Some("unrecognized arguments: --HELP")),
            (
                &["--help\u{e9}"],
                Some("unrecognized arguments: --help\u{e9}"),
            ),
        ]);
    }

    #[test]
    fn one_dash_before_an_equals_sign_begins_every_option_string() {
        assert_outcomes(&[
            (&["-=x"], Some("ambiguous option: -=x could match -h, --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size")),
            (&["-="], Some("ambiguous option: -= could match -h, --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size")),
            (&["-=x", "--assembly", "bad"], Some("ambiguous option: -=x could match -h, --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size")),
        ]);
    }

    #[test]
    fn of_two_defects_the_one_reached_first_is_reported() {
        assert_outcomes(&[
            (
                &["--assembly", "bad", "-=x"],
                Some("argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')"),
            ),
            (
                &["--assembly", "bad", "--r=1"],
                Some("argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')"),
            ),
            (
                &["--reference-max-sv-size", "abc", "--re"],
                Some("argument --reference-max-sv-size: invalid int value: 'abc'"),
            ),
        ]);
    }

    #[test]
    fn a_refused_value_is_quoted_as_python_writes_a_text() {
        assert_outcomes(&[
            (&["-h=it's"], Some("argument -h/--help: ignored explicit argument \"it's\"")),
            (&["--help=a\tb"], Some("argument -h/--help: ignored explicit argument 'a\\tb'")),
            (&["--help=\"'"], Some("argument -h/--help: ignored explicit argument '\"\\''")),
            (&["--assembly", "it's"], Some("argument --assembly: invalid choice: \"it's\" (choose from 'grch37', 'grch38')")),
            (&["--assembly=a\\b"], Some("argument --assembly: invalid choice: 'a\\\\b' (choose from 'grch37', 'grch38')")),
            (&["--reference-max-sv-size", "it's"], Some("argument --reference-max-sv-size: invalid int value: \"it's\"")),
            (&["--reference-max-sv-size=\u{e9}\n"], Some("argument --reference-max-sv-size: invalid int value: '\u{e9}\\n'")),
        ]);
    }

    #[test]
    fn a_byte_that_is_not_utf8_is_written_as_python_escapes_it() {
        let message = |arguments: &[&[u8]]| bytes(arguments).unwrap_err();
        assert_eq!(message(&[b"\xff"]), "unrecognized arguments: \\udcff");
        assert_eq!(
            message(&[b"x\xffy", b"-u\xc3"]),
            "unrecognized arguments: x\\udcffy -u\\udcc3"
        );
        assert_eq!(
            message(&[b"--", b"\xff"]),
            "unrecognized arguments: -- \\udcff"
        );
        assert_eq!(
            message(&[b"--re\xff"]),
            "unrecognized arguments: --re\\udcff"
        );
        assert_eq!(
            message(&[b"--assembly", b"\xff"]),
            "argument --assembly: invalid choice: '\\udcff' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            message(&[b"--assembly=gr\xe9"]),
            "argument --assembly: invalid choice: 'gr\\udce9' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            message(&[b"--re=\xff"]),
            "ambiguous option: --re=\\udcff could match --reference-release, --reference-max-sv-size"
        );
        assert_eq!(
            message(&[b"--help=\xff"]),
            "argument -h/--help: ignored explicit argument '\\udcff'"
        );
        assert_eq!(
            message(&[b"-h=\xe2\x82"]),
            "argument -h/--help: ignored explicit argument '\\udce2\\udc82'"
        );
        assert_eq!(
            message(&[b"--reference-max-sv-size", b"1\xff"]),
            "argument --reference-max-sv-size: invalid int value: '1\\udcff'"
        );
        assert_eq!(
            message(&[b"--perl-dir", b"-\xff"]),
            "argument --perl-dir: expected one argument"
        );
        assert_eq!(bytes(&[b"-h\xff"]), Ok(Request::Help));
        for value in [&b"-1\xff"[..], b"-\xd9\xa3\xff"] {
            let Ok(Request::Compare(options)) = bytes(&[b"--perl-dir", value]) else {
                panic!("a value that begins as a negative number: {value:?}");
            };
            assert_eq!(options.perl_dir, PathBuf::from(OsStr::from_bytes(value)));
        }
    }

    #[test]
    fn the_help_marks_the_default_release_and_no_other() {
        let marked: Vec<&str> = REFERENCE_RELEASES
            .into_iter()
            .filter(|release| HELP.contains(&format!("{release} (default)")))
            .collect();
        assert_eq!(marked, [DEFAULT_REFERENCE_RELEASE]);
        assert_eq!(HELP.matches("(default)").count(), 1);
    }

    #[test]
    fn the_help_is_one_text_that_opens_with_the_usage_block() {
        let lines: Vec<&str> = HELP.lines().collect();
        assert_eq!(lines.len(), 79);
        assert!(HELP.ends_with("route or dropped one that is not.\n"));
        assert!(lines.iter().all(|line| line.len() <= 78 && line.is_ascii()));
        assert_eq!(
            usage(),
            "usage: vep-sv-compare [-h] [--assembly {grch37,grch38}]\n\
             \x20                     [--input-dir INPUT_DIR] [--perl-dir PERL_DIR]\n\
             \x20                     [--rust-dir RUST_DIR] [--output-dir OUTPUT_DIR]\n\
             \x20                     [--vep-rs-cache VEP_RS_CACHE]\n\
             \x20                     [--scored-engine SCORED_ENGINE]\n\
             \x20                     [--reference-release {115.2,116.2,116.3}]\n\
             \x20                     [--reference-max-sv-size REFERENCE_MAX_SV_SIZE]\n"
        );
        assert_eq!(
            &lines[7..12],
            [
                "",
                "Compare Ensembl VEP and vep-rs outputs across a directory of SV VCFs.",
                "",
                "options:",
                "  -h, --help            show this help message and exit",
            ]
        );
        // Every option string heads its entry, in the order of the usage block.
        let entries: Vec<&str> = lines[11..]
            .iter()
            .filter(|line| line.starts_with("  -"))
            .map(|line| line.trim_start().split([' ', ',']).next().expect("a word"))
            .collect();
        assert_eq!(entries[0], OPTION_STRINGS[0]);
        assert_eq!(entries[1..], OPTION_STRINGS[2..]);
    }

    /// fastVEP writes a raw symbolic token in the Allele column through v0.5.0 and the
    /// class string itself from v0.6.0, so what the help says the mapping changes is
    /// bounded by those two versions.
    #[test]
    fn the_scored_engine_help_names_the_versions_whose_alleles_the_mapping_changes() {
        let text = HELP.split_whitespace().collect::<Vec<_>>().join(" ");
        for phrase in [
            "fastVEP through v0.5.0 writes the raw token, and without the mapping every",
            "tuple of such an output is a guaranteed miss; from v0.6.0 it writes the class \
             string itself, which the mapping leaves as written.",
            "the SO class string Ensembl VEP and vep-rs write there",
            "The reference side is never mapped.",
        ] {
            assert!(text.contains(phrase), "the help lacks {phrase:?}");
        }
        for phrase in [
            "tuple on that arm is a guaranteed miss",
            "the scored side's raw symbolic ALT tokens",
        ] {
            assert!(!text.contains(phrase), "the help holds {phrase:?}");
        }
    }

    /// The help calls the reference engine `Ensembl VEP` and a fault of vep-rs a defect.
    /// `perl` stands in it only where an option, that option's value or its default
    /// directory is named, each time all in small letters or all in capitals.
    #[test]
    fn the_help_names_the_reference_engine_ensembl_vep_and_a_vep_rs_fault_a_defect() {
        let text = HELP.split_whitespace().collect::<Vec<_>>().join(" ");
        for phrase in [
            "Compare Ensembl VEP and vep-rs outputs across a directory of SV VCFs.",
            "--perl-dir PERL_DIR Directory with Ensembl VEP output .txt files (default: \
             tmp/sv_perl/{assembly})",
            "chromosome mask, which removes Ensembl VEP tuples naming a transcript the \
             variant's chromosome does not contain.",
            "vep-rs naming a transcript the cache does not hold for the variant's chromosome \
             is a vep-rs defect and raises,",
        ] {
            assert!(text.contains(phrase), "the help lacks {phrase:?}");
        }
        assert!(
            !HELP.contains("Perl"),
            "the help names the engine by its language"
        );
        let refused: Vec<&str> = text
            .split(|c: char| !c.is_ascii_alphanumeric())
            .filter(|word| {
                ["bug", "bugs"]
                    .iter()
                    .any(|name| word.eq_ignore_ascii_case(name))
            })
            .collect();
        assert!(refused.is_empty(), "the help holds {refused:?}");
    }

    #[test]
    fn every_argument_from_the_first_double_dash_on_is_positional() {
        let message = |arguments: &[&str]| parse(arguments).unwrap_err();
        assert_eq!(message(&["--"]), "unrecognized arguments: --");
        assert_eq!(message(&["x", "--", "y"]), "unrecognized arguments: x -- y");
        assert_eq!(
            message(&["--", "--assembly", "bad"]),
            "unrecognized arguments: -- --assembly bad"
        );
        assert_eq!(
            message(&["--input-dir", "d", "--", "--assembly"]),
            "unrecognized arguments: -- --assembly"
        );
        assert_eq!(
            message(&["--", "--re", "1"]),
            "unrecognized arguments: -- --re 1"
        );
        assert_eq!(
            message(&["x", "--", "--re"]),
            "unrecognized arguments: x -- --re"
        );
        assert_eq!(
            message(&["--", "--output-dir", "out", "--", "--"]),
            "unrecognized arguments: -- --output-dir out -- --"
        );
        assert_eq!(
            message(&["--re", "1", "--", "x"]),
            "ambiguous option: --re could match --reference-release, --reference-max-sv-size"
        );
        assert_eq!(
            message(&["--assembly", "bad", "--", "x"]),
            "argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')"
        );
        assert_eq!(
            parse(&["--perl-dir=--"]).unwrap().perl_dir,
            PathBuf::from("--")
        );
    }

    #[test]
    fn two_dashes_before_an_equals_sign_begin_every_option() {
        assert_eq!(
            parse(&["--=x"]).unwrap_err(),
            "ambiguous option: --=x could match --help, --assembly, --input-dir, --perl-dir, --rust-dir, --output-dir, --vep-rs-cache, --scored-engine, --reference-release, --reference-max-sv-size"
        );
    }

    #[test]
    fn a_value_may_begin_with_a_dash_where_it_reads_as_no_option() {
        assert_eq!(
            parse(&["--perl-dir", "-1"]).unwrap().perl_dir,
            PathBuf::from("-1")
        );
        assert_eq!(
            parse(&["--perl-dir", "-.5"]).unwrap().perl_dir,
            PathBuf::from("-.5")
        );
        assert_eq!(
            parse(&["--perl-dir", "-x y"]).unwrap().perl_dir,
            PathBuf::from("-x y")
        );
        assert_eq!(
            parse(&["--perl-dir", "-"]).unwrap().perl_dir,
            PathBuf::from("-")
        );
        assert_eq!(
            parse(&["--perl-dir=-x"]).unwrap().perl_dir,
            PathBuf::from("-x")
        );
        assert_eq!(
            parse(&["--perl-dir", "-1e5"]).unwrap().perl_dir,
            PathBuf::from("-1e5")
        );
        assert_eq!(
            parse(&["--perl-dir", "--foo bar"]).unwrap().perl_dir,
            PathBuf::from("--foo bar")
        );
    }
}
