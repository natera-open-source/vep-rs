// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as Python's `argparse` reads two positional arguments, two
//! one-value options and `-h`: an option by its full name or any unambiguous prefix of it,
//! its value as the next argument or after `=`, a value that starts with `-` taken only
//! where it reads as a negative number or carries a space, the last of a repeated option
//! kept, the positional arguments before, between or after the options, every argument
//! after the first `--` a value, help where `-h` or `--help` is met before the first fault.

use crate::pyfmt::{py_repr, shown, surrogateescape};
use crate::unicode::decimal_value;

pub const USAGE: &str =
    "usage: vep-fastvep-plugin-tuples [-h] [--field-map FIELD_MAP] [--fields FIELDS] input tuples";
pub const HELP: &str = "\
Write the supplementary columns of fastVEP's tab output as plugin tuples.

positional arguments:
  input                 fastVEP tab output, as written or as normalized
  tuples                tuple file to write

options:
  -h, --help            show this help message and exit
  --field-map FIELD_MAP
                        field map (default: data/fastvep_field_map.tsv of this tool)
  --fields FIELDS       comma-separated Ensembl VEP field names to keep";

/// The arguments of one run, each the bytes of its argument as written.
#[derive(Debug, PartialEq, Eq)]
pub struct Args {
    pub input: Vec<u8>,
    pub tuples: Vec<u8>,
    pub field_map: Option<Vec<u8>>,
    pub fields: Option<Vec<u8>>,
}

/// What the command line asks for.
#[derive(Debug, PartialEq, Eq)]
pub enum Parsed {
    Run(Args),
    Help,
    /// A usage error, with its message.
    Usage(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Action {
    Help,
    FieldMap,
    Fields,
}

/// The option strings in the order the parser declares them.
const OPTIONS: [(&str, Action); 4] = [
    ("-h", Action::Help),
    ("--help", Action::Help),
    ("--field-map", Action::FieldMap),
    ("--fields", Action::Fields),
];

/// The positional arguments in the order the parser declares them.
const POSITIONALS: [&str; 2] = ["input", "tuples"];

fn action_of(option: &[u8]) -> Option<Action> {
    OPTIONS
        .iter()
        .find(|(name, _)| name.as_bytes() == option)
        .map(|(_, action)| *action)
}

/// What joins an option to a value written in the same argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Joint {
    /// `--option=value`.
    Equals,
    /// `-hvalue`.
    Nothing,
}

/// One reading of an argument as an option: the action (None for an option this parser
/// does not have), the option string and the value written in the same argument.
#[derive(Clone, Debug)]
struct Reading {
    action: Option<Action>,
    option: Vec<u8>,
    explicit: Option<(Joint, Vec<u8>)>,
}

fn split_at_equals(arg: &[u8]) -> (&[u8], Option<&[u8]>) {
    match arg.iter().position(|b| *b == b'=') {
        Some(at) => (&arg[..at], Some(&arg[at + 1..])),
        None => (arg, None),
    }
}

/// Every option an argument that is no option string itself may stand for: a long option
/// by a prefix of its name, `-h` by the argument's first two characters.
fn readings_by_prefix(arg: &[u8]) -> Vec<Reading> {
    let (prefix, value) = split_at_equals(arg);
    let explicit = value.map(|value| (Joint::Equals, value.to_vec()));
    let mut readings = Vec::new();
    for (name, action) in OPTIONS {
        let name = name.as_bytes();
        if arg[1] != b'-' && name == &arg[..2] {
            readings.push(Reading {
                action: Some(action),
                option: name.to_vec(),
                explicit: Some((Joint::Nothing, arg[2..].to_vec())),
            });
        } else if name.starts_with(prefix) {
            readings.push(Reading {
                action: Some(action),
                option: name.to_vec(),
                explicit: explicit.clone(),
            });
        }
    }
    readings
}

/// Whether the argument opens like a negative number: `-`, an optional `.`, a decimal
/// digit of any script.
fn opens_like_a_negative_number(arg: &[u8]) -> bool {
    let rest = arg.strip_prefix(b"-").unwrap_or(arg);
    let rest = rest.strip_prefix(b".").unwrap_or(rest);
    let first = String::from_utf8_lossy(&rest[..rest.len().min(4)])
        .chars()
        .next();
    first.is_some_and(|c| decimal_value(c).is_some())
}

/// The readings of an argument as an option, or None where it is a value.
fn readings(arg: &[u8]) -> Option<Vec<Reading>> {
    if arg.first() != Some(&b'-') {
        return None;
    }
    let exact = |option: &[u8], explicit: Option<(Joint, Vec<u8>)>| {
        action_of(option).map(|action| {
            vec![Reading {
                action: Some(action),
                option: option.to_vec(),
                explicit,
            }]
        })
    };
    if let Some(found) = exact(arg, None) {
        return Some(found);
    }
    if arg.len() == 1 {
        return None;
    }
    if let (option, Some(value)) = split_at_equals(arg) {
        if let Some(found) = exact(option, Some((Joint::Equals, value.to_vec()))) {
            return Some(found);
        }
    }
    let by_prefix = readings_by_prefix(arg);
    if !by_prefix.is_empty() {
        return Some(by_prefix);
    }
    if opens_like_a_negative_number(arg) || arg.contains(&b' ') {
        return None;
    }
    Some(vec![Reading {
        action: None,
        option: arg.to_vec(),
        explicit: None,
    }])
}

/// One argument's place in the command line.
enum Kind {
    Value,
    /// The first `--`; every argument after it is a value.
    Separator,
    Option(Vec<Reading>),
}

/// The positional arguments still to be read and the values of those that are.
struct Positionals<'a> {
    values: Vec<&'a [u8]>,
}

impl<'a> Positionals<'a> {
    /// Reads as many of the positional arguments still open as the values from `start` on
    /// fill, each one value with the `--` beside it, and returns where the reading stopped.
    fn consume(&mut self, argv: &'a [Vec<u8>], kinds: &[Kind], start: usize) -> usize {
        let open = POSITIONALS.len() - self.values.len();
        // The lengths of `(-*A-*)` matched `count` times from `start`, or None.
        let matched = |count: usize| {
            let mut at = start;
            let mut spans = Vec::with_capacity(count);
            for _ in 0..count {
                let from = at;
                while matches!(kinds.get(at), Some(Kind::Separator)) {
                    at += 1;
                }
                if !matches!(kinds.get(at), Some(Kind::Value)) {
                    return None;
                }
                let value = at;
                at += 1;
                while matches!(kinds.get(at), Some(Kind::Separator)) {
                    at += 1;
                }
                spans.push((from, value, at));
            }
            Some(spans)
        };
        let Some(spans) = (1..=open).rev().find_map(matched) else {
            return start;
        };
        let mut end = start;
        for (_, value, after) in spans {
            self.values.push(&argv[value]);
            end = after;
        }
        end
    }
}

/// Reads the arguments after the program's name.
pub fn parse(argv: &[Vec<u8>]) -> Parsed {
    let mut kinds = Vec::with_capacity(argv.len());
    let mut separated = false;
    for arg in argv {
        kinds.push(if separated {
            Kind::Value
        } else if arg == b"--" {
            separated = true;
            Kind::Separator
        } else {
            readings(arg).map_or(Kind::Value, Kind::Option)
        });
    }
    let is_option = |at: usize| matches!(kinds.get(at), Some(Kind::Option(_)));
    let last_option = (0..argv.len()).rev().find(|at| is_option(*at));
    let mut positionals = Positionals { values: Vec::new() };
    let mut field_map = None;
    let mut fields = None;
    // A value no positional argument and no option takes is unrecognized; the fault is
    // raised after every argument is read, and after a missing positional argument.
    let mut unrecognized: Vec<String> = Vec::new();
    let mut at = 0;
    while last_option.is_some_and(|last| at <= last) {
        let next_option = (at..argv.len())
            .find(|index| is_option(*index))
            .expect("an option at or after this argument");
        if at != next_option {
            let end = positionals.consume(argv, &kinds, at);
            if end > at {
                at = end;
                continue;
            }
            unrecognized.extend(argv[at..next_option].iter().map(|arg| shown(arg)));
            at = next_option;
        }
        let Kind::Option(readings) = &kinds[at] else {
            unreachable!("the next option is an option");
        };
        if readings.len() > 1 {
            let matches: Vec<String> = readings.iter().map(|r| shown(&r.option)).collect();
            return Parsed::Usage(format!(
                "ambiguous option: {} could match {}",
                shown(&argv[at]),
                matches.join(", ")
            ));
        }
        let first = &readings[0];
        let mut option: &[u8] = &first.option;
        let mut explicit: Option<(Joint, &[u8])> = first
            .explicit
            .as_ref()
            .map(|(joint, value)| (*joint, value.as_slice()));
        let mut help = false;
        let mut stored: Option<(Action, &[u8])> = None;
        let mut next = at + 1;
        let ignored = |value: &[u8]| {
            Parsed::Usage(format!(
                "argument -h/--help: ignored explicit argument {}",
                py_repr(&surrogateescape(value))
            ))
        };
        match first.action {
            None => unrecognized.push(shown(&argv[at])),
            Some(Action::Help) => loop {
                match explicit {
                    None => {
                        help = true;
                        break;
                    }
                    // `-h` run together with more: each further character is a flag of
                    // its own, and no other flag exists.
                    Some((joint, tail)) if option == b"-h" && !tail.is_empty() => {
                        if joint == Joint::Equals || tail[0] == b'-' {
                            return ignored(tail);
                        }
                        help = true;
                        if tail[0] != b'h' {
                            unrecognized.push(format!("-{}", shown(tail)));
                            break;
                        }
                        option = b"-h";
                        explicit = match &tail[1..] {
                            [] => None,
                            [b'=', value @ ..] => Some((Joint::Equals, value)),
                            value => Some((Joint::Nothing, value)),
                        };
                    }
                    Some((_, value)) => return ignored(value),
                }
            },
            Some(action) => match explicit {
                Some((_, value)) => stored = Some((action, value)),
                None if matches!(kinds.get(next), Some(Kind::Value)) => {
                    stored = Some((action, &argv[next]));
                    next += 1;
                }
                None => {
                    return Parsed::Usage(format!(
                        "argument {}: expected one argument",
                        shown(option)
                    ));
                }
            },
        }
        if help {
            return Parsed::Help;
        }
        match stored {
            Some((Action::FieldMap, value)) => field_map = Some(value.to_vec()),
            Some((Action::Fields, value)) => fields = Some(value.to_vec()),
            Some((Action::Help, _)) => unreachable!("help stores no value"),
            None => {}
        }
        at = next;
    }
    let stop = positionals.consume(argv, &kinds, at);
    unrecognized.extend(argv[stop..].iter().map(|arg| shown(arg)));
    if positionals.values.len() < POSITIONALS.len() {
        return Parsed::Usage(format!(
            "the following arguments are required: {}",
            POSITIONALS[positionals.values.len()..].join(", ")
        ));
    }
    if !unrecognized.is_empty() {
        return Parsed::Usage(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }
    Parsed::Run(Args {
        input: positionals.values[0].to_vec(),
        tuples: positionals.values[1].to_vec(),
        field_map,
        fields,
    })
}

#[cfg(test)]
#[path = "cli_argparse_tests.rs"]
mod argparse_tests;
