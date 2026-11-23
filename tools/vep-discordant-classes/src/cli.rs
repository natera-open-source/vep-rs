// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as CPython's `argparse` reads a parser of one positional
//! (`input`), `-o/--output` and `-h/--help` with its defaults. Every argument is first read
//! on its own as a value, as the `--` marker or as an option; the arguments are then taken
//! from left to right, the input by the first value that no option takes. Errors come in
//! argparse's order: what an option meets as it is read, then the missing input, then the
//! arguments nothing took.

use crate::pyfmt::is_decimal_digit;
use memchr::memchr;
use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;

pub const USAGE: &str = "usage: vep-discordant-classes [-h] [-o OUTPUT] input";
pub const HELP: &str = "\
Classify discordant tuples

positional arguments:
  input                Path to discordant.tsv

options:
  -h, --help           show this help message and exit
  -o, --output OUTPUT  Output path (default: discordant_classified.tsv)";

#[derive(Debug, PartialEq, Eq)]
pub struct Args {
    pub input: OsString,
    /// None where the option is not given: the table is then written beside the input.
    pub output: Option<OsString>,
}

/// What a command line asks for.
#[derive(Debug, PartialEq, Eq)]
pub enum Command {
    Run(Args),
    Help,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Action {
    Help,
    Output,
}

/// The option strings in the order an ambiguity lists them.
const OPTION_STRINGS: [(&str, Action); 4] = [
    ("-h", Action::Help),
    ("--help", Action::Help),
    ("-o", Action::Output),
    ("--output", Action::Output),
];

/// One way an argument reads as an option: the option string it names and the value it
/// carries in the same argument, with whether a `=` stands before that value.
#[derive(Clone, Copy)]
struct Reading<'a> {
    action: Action,
    option: &'static str,
    explicit: Option<(bool, &'a [u8])>,
}

/// How one argument reads before any is consumed.
enum Kind<'a> {
    /// A value: the input, or what `-o` takes.
    Argument,
    /// The first `--`: every argument after it is a value.
    Marker,
    /// Every option the argument can name; more than one is an ambiguity.
    Options(Vec<Reading<'a>>),
    /// Shaped like an option and naming none.
    Unknown,
}

fn lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Whether the argument starts as the pattern `-\.?\d` does: such an argument is a value
/// where it names no option.
fn starts_like_a_negative_number(arg: &[u8]) -> bool {
    let rest = &arg[1..];
    let rest = rest.strip_prefix(b".").unwrap_or(rest);
    (1..=rest.len().min(4))
        .find_map(|length| std::str::from_utf8(&rest[..length]).ok())
        .and_then(|text| text.chars().next())
        .is_some_and(is_decimal_digit)
}

/// Reads one argument that stands before any `--` marker.
fn read(arg: &[u8]) -> Kind<'_> {
    if arg.first() != Some(&b'-') {
        return Kind::Argument;
    }
    let named = |name: &[u8]| {
        OPTION_STRINGS
            .iter()
            .find(|(option, _)| option.as_bytes() == name)
    };
    if let Some((option, action)) = named(arg) {
        return Kind::Options(vec![Reading {
            action: *action,
            option,
            explicit: None,
        }]);
    }
    if arg.len() == 1 {
        return Kind::Argument;
    }
    let (prefix, after_equals) = match memchr(b'=', arg) {
        Some(at) => (&arg[..at], Some((true, &arg[at + 1..]))),
        None => (arg, None),
    };
    if let (Some((option, action)), Some(_)) = (named(prefix), after_equals) {
        return Kind::Options(vec![Reading {
            action: *action,
            option,
            explicit: after_equals,
        }]);
    }
    let mut readings = Vec::new();
    for (option, action) in OPTION_STRINGS {
        // A one-letter option takes the rest of its argument as its value; any option is
        // named by a prefix that it starts with.
        let explicit = if arg[1] != b'-' && option.as_bytes() == &arg[..2] {
            Some((false, &arg[2..]))
        } else if option.as_bytes().starts_with(prefix) {
            after_equals
        } else {
            continue;
        };
        readings.push(Reading {
            action,
            option,
            explicit,
        });
    }
    if !readings.is_empty() {
        Kind::Options(readings)
    } else if starts_like_a_negative_number(arg) || arg.contains(&b' ') {
        Kind::Argument
    } else {
        Kind::Unknown
    }
}

/// Takes the input from the arguments at `at` unless one is taken already: one value, with
/// the marker before or after it. Returns the index behind what it took.
fn take_input<'a>(
    at: usize,
    kinds: &[Kind<'_>],
    argv: &'a [OsString],
    input: &mut Option<&'a OsStr>,
) -> usize {
    if input.is_some() {
        return at;
    }
    let mut next = at;
    if matches!(kinds.get(next), Some(Kind::Marker)) {
        next += 1;
    }
    if !matches!(kinds.get(next), Some(Kind::Argument)) {
        return at;
    }
    *input = Some(&argv[next]);
    next += 1;
    if matches!(kinds.get(next), Some(Kind::Marker)) {
        next += 1;
    }
    next
}

/// What the option at `at` does: the index behind what it took, or None where it asks for
/// the help. A `-h` reads the characters behind it as further one-letter options.
fn take_option<'a>(
    at: usize,
    reading: Reading<'a>,
    kinds: &[Kind<'_>],
    argv: &'a [OsString],
    output: &mut Option<&'a OsStr>,
    unread: &mut Vec<String>,
) -> Result<Option<usize>, String> {
    let Reading {
        mut action,
        mut option,
        mut explicit,
    } = reading;
    let mut help = false;
    let stop = loop {
        match (action, explicit) {
            (Action::Output, Some((_, value))) => {
                *output = Some(OsStr::from_bytes(value));
                break at + 1;
            }
            (Action::Output, None) => {
                if !matches!(kinds.get(at + 1), Some(Kind::Argument)) {
                    return Err("argument -o/--output: expected one argument".to_string());
                }
                *output = Some(&argv[at + 1]);
                break at + 2;
            }
            (Action::Help, None) => {
                help = true;
                break at + 1;
            }
            (Action::Help, Some((after_equals, value))) => {
                let refused = after_equals || value.first() == Some(&b'-');
                if option != "-h" || value.is_empty() || refused {
                    return Err(format!(
                        "argument -h/--help: ignored explicit argument '{}'",
                        lossy(value)
                    ));
                }
                help = true;
                (action, option) = match value[0] {
                    b'h' => (Action::Help, "-h"),
                    b'o' => (Action::Output, "-o"),
                    _ => {
                        unread.push(format!("-{}", lossy(value)));
                        break at + 1;
                    }
                };
                explicit = match &value[1..] {
                    [] => None,
                    [b'=', rest @ ..] => Some((true, rest)),
                    rest => Some((false, rest)),
                };
            }
        }
    };
    Ok(if help { None } else { Some(stop) })
}

pub fn parse(argv: &[OsString]) -> Result<Command, String> {
    let mut kinds: Vec<Kind<'_>> = Vec::with_capacity(argv.len());
    let mut after_marker = false;
    for arg in argv {
        let arg = arg.as_bytes();
        kinds.push(if after_marker {
            Kind::Argument
        } else if arg == b"--" {
            after_marker = true;
            Kind::Marker
        } else {
            read(arg)
        });
    }
    let is_option = |index: &usize| matches!(kinds[*index], Kind::Options(_) | Kind::Unknown);
    let mut input = None;
    let mut output = None;
    let mut unread: Vec<String> = Vec::new();
    let mut at = 0;
    if let Some(last) = (0..argv.len()).rev().find(is_option) {
        while at <= last {
            let next = (at..=last).find(is_option).unwrap_or(last);
            if at != next {
                let behind = take_input(at, &kinds, argv, &mut input);
                if behind > at {
                    at = behind;
                    continue;
                }
                unread.extend(argv[at..next].iter().map(|arg| lossy(arg.as_bytes())));
                at = next;
            }
            let Kind::Options(readings) = &kinds[at] else {
                unread.push(lossy(argv[at].as_bytes()));
                at += 1;
                continue;
            };
            if let [_, _, ..] = readings.as_slice() {
                let matches: Vec<&str> = readings.iter().map(|reading| reading.option).collect();
                return Err(format!(
                    "ambiguous option: {} could match {}",
                    lossy(argv[at].as_bytes()),
                    matches.join(", ")
                ));
            }
            match take_option(at, readings[0], &kinds, argv, &mut output, &mut unread)? {
                Some(behind) => at = behind,
                None => return Ok(Command::Help),
            }
        }
    }
    let behind = take_input(at, &kinds, argv, &mut input);
    unread.extend(argv[behind..].iter().map(|arg| lossy(arg.as_bytes())));
    match input {
        None => Err("the following arguments are required: input".to_string()),
        Some(_) if !unread.is_empty() => {
            Err(format!("unrecognized arguments: {}", unread.join(" ")))
        }
        Some(input) => Ok(Command::Run(Args {
            input: input.to_os_string(),
            output: output.map(OsStr::to_os_string),
        })),
    }
}

#[cfg(test)]
#[path = "cli_parse_tests.rs"]
mod parse_tests;
