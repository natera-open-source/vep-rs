// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A command line of long options that each take one value, read as CPython's
//! `argparse` reads it: an option may be shortened to any prefix that fits it alone,
//! takes its value as the next argument or after `=`, and a later value replaces an
//! earlier one. A value that starts with `-` is taken only when it is `-` alone, starts
//! like a negative number or holds a space. `-h` and `--help` ask for the usage.

use std::collections::HashMap;

/// One option of the command line.
#[derive(Debug, Clone, Copy)]
pub struct Opt {
    /// The option with its two dashes.
    pub name: &'static str,
    pub required: bool,
}

/// What the command line asks for.
#[derive(Debug, PartialEq, Eq)]
pub enum Parsed {
    /// The value of every option given, by its name.
    Values(HashMap<&'static str, Vec<u8>>),
    Help,
}

/// How one argument reads before any is consumed.
enum Reading<'a> {
    /// A value for the option before it, or a stray positional.
    Argument,
    /// One option, with the text after `=` when the argument carries one.
    Option(&'static str, Option<&'a [u8]>),
    /// `-h` with more characters behind it.
    ShortHelp(&'a [u8]),
    /// A prefix that more than one option starts with.
    Ambiguous(Vec<&'static str>),
    /// Shaped like an option and naming none.
    Unknown,
}

const HELP: &str = "--help";

fn read<'a>(options: &[Opt], arg: &'a [u8]) -> Reading<'a> {
    if arg.first() != Some(&b'-') || arg.len() == 1 {
        return Reading::Argument;
    }
    let (head, inline) = match arg.iter().position(|byte| *byte == b'=') {
        Some(at) => (&arg[..at], Some(&arg[at + 1..])),
        None => (arg, None),
    };
    if head == b"-h" {
        return Reading::Option(HELP, inline);
    }
    if arg[1] == b'-' {
        let names = || options.iter().map(|option| option.name).chain([HELP]);
        let matches: Vec<&'static str> = match names().find(|name| name.as_bytes() == head) {
            Some(name) => vec![name],
            None => names()
                .filter(|name| name.as_bytes().starts_with(head))
                .collect(),
        };
        match matches.as_slice() {
            [] => {}
            [name] => return Reading::Option(name, inline),
            _ => return Reading::Ambiguous(matches),
        }
    } else if arg.starts_with(b"-h") {
        return Reading::ShortHelp(&arg[2..]);
    }
    let negative_number = matches!(arg, [b'-', b'0'..=b'9', ..] | [b'-', b'.', b'0'..=b'9', ..]);
    if negative_number || arg.contains(&b' ') {
        Reading::Argument
    } else {
        Reading::Unknown
    }
}

fn lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The characters behind a `-h`: a further `h` is help again, a letter that names no
/// option is left unread and the help is shown all the same, and a `-` or an `=` there
/// is a value the option refuses.
fn short_help(mut rest: &[u8]) -> Result<(), String> {
    let refused = |value: &[u8]| {
        format!(
            "argument -h/--help: ignored explicit argument '{}'",
            lossy(value)
        )
    };
    loop {
        match rest {
            [b'-', ..] => return Err(refused(rest)),
            [b'h', b'=', value @ ..] => return Err(refused(value)),
            [b'h', more @ ..] if !more.is_empty() => rest = more,
            _ => return Ok(()),
        }
    }
}

/// Reads the arguments after the program's name. The error is argparse's message.
pub fn parse(options: &[Opt], argv: &[Vec<u8>]) -> Result<Parsed, String> {
    let mut values: HashMap<&'static str, Vec<u8>> = HashMap::new();
    let mut unread: Vec<String> = Vec::new();
    let mut index = 0;
    while index < argv.len() {
        let arg = argv[index].as_slice();
        index += 1;
        if arg == b"--" {
            // The marker makes every later argument a positional, and the parser takes none.
            unread.extend(argv[index - 1..].iter().map(|arg| lossy(arg)));
            break;
        }
        let (name, inline) = match read(options, arg) {
            Reading::Argument | Reading::Unknown => {
                unread.push(lossy(arg));
                continue;
            }
            Reading::Ambiguous(matches) => {
                return Err(format!(
                    "ambiguous option: {} could match {}",
                    lossy(arg),
                    matches.join(", ")
                ))
            }
            Reading::ShortHelp(rest) => return short_help(rest).map(|()| Parsed::Help),
            Reading::Option(name, inline) => (name, inline),
        };
        if name == HELP {
            return match inline {
                Some(value) => Err(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    lossy(value)
                )),
                None => Ok(Parsed::Help),
            };
        }
        let value = match inline {
            Some(value) => value,
            None => match argv.get(index).map(Vec::as_slice) {
                Some(next) if next != b"--" && matches!(read(options, next), Reading::Argument) => {
                    index += 1;
                    next
                }
                _ => return Err(format!("argument {name}: expected one argument")),
            },
        };
        values.insert(name, value.to_vec());
    }
    let missing: Vec<&str> = options
        .iter()
        .filter(|option| option.required && !values.contains_key(option.name))
        .map(|option| option.name)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    if !unread.is_empty() {
        return Err(format!("unrecognized arguments: {}", unread.join(" ")));
    }
    Ok(Parsed::Values(values))
}

#[cfg(test)]
mod tests {
    use super::{parse, Opt, Parsed};

    const OPTIONS: [Opt; 4] = [
        Opt {
            name: "--version-a",
            required: true,
        },
        Opt {
            name: "--version-b",
            required: true,
        },
        Opt {
            name: "--out-dir",
            required: false,
        },
        Opt {
            name: "--sweep-id",
            required: false,
        },
    ];

    fn run(arguments: &[&str]) -> Result<Parsed, String> {
        let argv: Vec<Vec<u8>> = arguments
            .iter()
            .map(|argument| argument.as_bytes().to_vec())
            .collect();
        parse(&OPTIONS, &argv)
    }

    fn value(parsed: &Result<Parsed, String>, name: &str) -> Option<String> {
        match parsed {
            Ok(Parsed::Values(values)) => values
                .get(name)
                .map(|value| String::from_utf8_lossy(value).into_owned()),
            _ => None,
        }
    }

    /// Each expectation is what `argparse.ArgumentParser` of CPython 3.14 does with the
    /// same four options and arguments.
    #[test]
    fn values_follow_their_option_or_an_equals_sign() {
        let parsed = run(&["--version-a", "113", "--version-b=115", "--out-dir", "o"]);
        assert_eq!(value(&parsed, "--version-a").as_deref(), Some("113"));
        assert_eq!(value(&parsed, "--version-b").as_deref(), Some("115"));
        assert_eq!(value(&parsed, "--out-dir").as_deref(), Some("o"));
        assert_eq!(value(&parsed, "--sweep-id"), None);
        // A later value replaces an earlier one; an empty value after `=` is a value.
        let parsed = run(&["--version-a=1", "--version-b", "2", "--version-a="]);
        assert_eq!(value(&parsed, "--version-a").as_deref(), Some(""));
    }

    #[test]
    fn a_prefix_that_fits_one_option_names_it() {
        let parsed = run(&[
            "--version-a",
            "1",
            "--version-b",
            "2",
            "--o",
            "dir",
            "--s=x",
        ]);
        assert_eq!(value(&parsed, "--out-dir").as_deref(), Some("dir"));
        assert_eq!(value(&parsed, "--sweep-id").as_deref(), Some("x"));
        assert_eq!(
            run(&["--version", "1"]),
            Err("ambiguous option: --version could match --version-a, --version-b".to_string())
        );
    }

    #[test]
    fn a_value_that_starts_with_a_dash_is_taken_when_it_is_no_option() {
        let parsed = run(&["--version-a", "-1", "--version-b", "-", "--out-dir", "-a b"]);
        assert_eq!(value(&parsed, "--version-a").as_deref(), Some("-1"));
        assert_eq!(value(&parsed, "--version-b").as_deref(), Some("-"));
        assert_eq!(value(&parsed, "--out-dir").as_deref(), Some("-a b"));
        assert_eq!(
            run(&["--version-a", "-x", "--version-b", "2"]),
            Err("argument --version-a: expected one argument".to_string())
        );
        assert_eq!(
            run(&["--version-a", "1", "--version-b"]),
            Err("argument --version-b: expected one argument".to_string())
        );
    }

    #[test]
    fn what_is_missing_is_named_before_what_is_not_read() {
        assert_eq!(
            run(&["stray"]),
            Err("the following arguments are required: --version-a, --version-b".to_string())
        );
        assert_eq!(
            run(&["--version-a", "1", "--version-b", "2", "stray", "--unknown"]),
            Err("unrecognized arguments: stray --unknown".to_string())
        );
        assert_eq!(
            run(&["--version-a", "1", "--version-b", "2", "--", "x"]),
            Err("unrecognized arguments: -- x".to_string())
        );
    }

    #[test]
    fn help_is_asked_for_in_several_ways() {
        for arguments in [
            &["-h"][..],
            &["--help"],
            &["--he"],
            &["-hh"],
            &["--version-a", "1", "-h"],
        ] {
            assert_eq!(run(arguments), Ok(Parsed::Help), "{arguments:?}");
        }
        assert_eq!(
            run(&["--help=x"]),
            Err("argument -h/--help: ignored explicit argument 'x'".to_string())
        );
    }
}
