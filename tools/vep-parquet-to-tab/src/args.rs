// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as Python's `argparse` reads it: one optional positional path,
//! `--input`/`-i`, `--output`/`-o`, `--canonical`, `--types-of` and `--help`/`-h`; a long
//! option under any prefix that names it alone, `--name=value`, `-ivalue` and `-i=value`;
//! `--` before arguments that are all positional; the last value of a repeated option.

use std::ffi::OsString;
use std::os::unix::ffi::{OsStrExt, OsStringExt};

/// What the command line asks for.
#[derive(Debug, PartialEq, Eq)]
pub struct Args {
    /// The output to read: `--input` when it is given and not empty, else the positional.
    pub path: Vec<u8>,
    /// The file to write, `-` for standard output.
    pub output: Vec<u8>,
    pub canonical: bool,
    /// The Parquet output whose stored types the canonical form of a `--tab` output
    /// follows, when one is named and not empty.
    pub types_of: Option<Vec<u8>>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Parsed {
    Run(Args),
    Help,
    /// A usage error: the message after `error: `.
    Usage(String),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Flag {
    Help,
    Input,
    Output,
    Canonical,
    TypesOf,
}

const LONG: [(&[u8], Flag); 5] = [
    (b"--help", Flag::Help),
    (b"--input", Flag::Input),
    (b"--output", Flag::Output),
    (b"--canonical", Flag::Canonical),
    (b"--types-of", Flag::TypesOf),
];
const SHORT: [(&[u8], Flag); 3] = [
    (b"-h", Flag::Help),
    (b"-i", Flag::Input),
    (b"-o", Flag::Output),
];

impl Flag {
    fn names(self) -> &'static str {
        match self {
            Flag::Help => "-h/--help",
            Flag::Input => "--input/-i",
            Flag::Output => "--output/-o",
            Flag::Canonical => "--canonical",
            Flag::TypesOf => "--types-of",
        }
    }
}

/// How one argument reads before any `--`.
enum Kind {
    Positional,
    /// An option, with the value written into the same argument. `attached` is a value
    /// that follows a short option with no `=` between them.
    Option {
        flag: Flag,
        value: Option<Vec<u8>>,
        attached: bool,
    },
    /// Looks like an option and names none.
    Unknown,
    Ambiguous(String),
}

fn exact(arg: &[u8]) -> Option<Flag> {
    LONG.iter()
        .chain(&SHORT)
        .find(|(name, _)| *name == arg)
        .map(|(_, flag)| *flag)
}

fn lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn classify(arg: &[u8]) -> Kind {
    if arg.is_empty() || arg[0] != b'-' {
        return Kind::Positional;
    }
    let option = |flag, value: Option<&[u8]>, attached| Kind::Option {
        flag,
        value: value.map(<[u8]>::to_vec),
        attached,
    };
    if let Some(flag) = exact(arg) {
        return option(flag, None, false);
    }
    if arg.len() == 1 {
        return Kind::Positional;
    }
    let (name, value) = match arg.iter().position(|b| *b == b'=') {
        Some(at) => (&arg[..at], Some(&arg[at + 1..])),
        None => (arg, None),
    };
    if let (Some(flag), Some(_)) = (exact(name), value) {
        return option(flag, value, false);
    }
    if arg[1] == b'-' {
        let matches: Vec<&(&[u8], Flag)> = LONG
            .iter()
            .filter(|(long, _)| long.starts_with(name))
            .collect();
        match matches.as_slice() {
            [] => {}
            [(_, flag)] => return option(*flag, value, false),
            several => {
                let names: Vec<String> = several.iter().map(|(long, _)| lossy(long)).collect();
                return Kind::Ambiguous(format!(
                    "ambiguous option: {} could match {}",
                    lossy(arg),
                    names.join(", ")
                ));
            }
        }
    } else if let Some((_, flag)) = SHORT.iter().find(|(short, _)| *short == &arg[..2]) {
        return option(*flag, Some(&arg[2..]), true);
    }
    // A negative number and an argument holding a space are positional.
    let numeric = match &arg[1..] {
        [b'.', digit, ..] | [digit, ..] => digit.is_ascii_digit(),
        [] => false,
    };
    if numeric || arg.contains(&b' ') {
        return Kind::Positional;
    }
    Kind::Unknown
}

/// Reads the arguments after the program name.
pub fn parse(argv: Vec<OsString>) -> Parsed {
    let argv: Vec<Vec<u8>> = argv.into_iter().map(OsStringExt::into_vec).collect();
    let separator = argv.iter().position(|arg| arg == b"--");
    let options_end = separator.unwrap_or(argv.len());
    let mut kinds = Vec::with_capacity(options_end);
    for arg in &argv[..options_end] {
        match classify(arg) {
            Kind::Ambiguous(message) => return Parsed::Usage(message),
            kind => kinds.push(kind),
        }
    }
    let mut positionals: Vec<&[u8]> = Vec::new();
    let mut unknown: Vec<&[u8]> = Vec::new();
    let mut input: Option<Vec<u8>> = None;
    let mut output = b"-".to_vec();
    let mut canonical = false;
    let mut types_of: Option<Vec<u8>> = None;
    let mut index = 0;
    while index < options_end {
        let arg = &argv[index];
        index += 1;
        match &kinds[index - 1] {
            Kind::Positional => positionals.push(arg),
            Kind::Unknown => unknown.push(arg),
            Kind::Ambiguous(_) => unreachable!("an ambiguous option ends the parse"),
            Kind::Option {
                flag,
                value,
                attached,
            } => match flag {
                Flag::Help | Flag::Canonical => {
                    if let Some(value) = value.as_ref().filter(|_| !attached) {
                        return Parsed::Usage(format!(
                            "argument {}: ignored explicit argument '{}'",
                            flag.names(),
                            lossy(value)
                        ));
                    }
                    if *flag == Flag::Help {
                        return Parsed::Help;
                    }
                    canonical = true;
                }
                Flag::Input | Flag::Output | Flag::TypesOf => {
                    let value = match value {
                        Some(value) => value.clone(),
                        None => {
                            let next = (index < options_end
                                && matches!(kinds[index], Kind::Positional))
                            .then(|| argv[index].clone());
                            let Some(next) = next else {
                                return Parsed::Usage(format!(
                                    "argument {}: expected one argument",
                                    flag.names()
                                ));
                            };
                            index += 1;
                            next
                        }
                    };
                    match flag {
                        Flag::Input => input = Some(value),
                        Flag::TypesOf => types_of = Some(value),
                        _ => output = value,
                    }
                }
            },
        }
    }
    if let Some(separator) = separator {
        positionals.extend(argv[separator + 1..].iter().map(Vec::as_slice));
    }
    unknown.extend(positionals.iter().skip(1));
    if !unknown.is_empty() {
        let unknown: Vec<String> = unknown.iter().map(|arg| lossy(arg)).collect();
        return Parsed::Usage(format!("unrecognized arguments: {}", unknown.join(" ")));
    }
    let path = input
        .filter(|path| !path.is_empty())
        .or_else(|| positionals.first().map(|path| path.to_vec()))
        .filter(|path| !path.is_empty());
    match path {
        Some(path) => Parsed::Run(Args {
            path,
            output,
            canonical,
            types_of: types_of.filter(|named| !named.is_empty()),
        }),
        None => Parsed::Usage("a Parquet path is required".to_string()),
    }
}

/// The argument's bytes as an `OsString`.
pub fn os(bytes: &[u8]) -> OsString {
    std::ffi::OsStr::from_bytes(bytes).to_os_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(args: &[&str]) -> Parsed {
        parse(args.iter().map(OsString::from).collect())
    }

    fn run(args: &[&str]) -> (String, String, bool) {
        match parsed(args) {
            Parsed::Run(args) => (lossy(&args.path), lossy(&args.output), args.canonical),
            other => panic!("{args:?}: {other:?}"),
        }
    }

    fn usage(args: &[&str]) -> String {
        match parsed(args) {
            Parsed::Usage(message) => message,
            other => panic!("{args:?}: {other:?}"),
        }
    }

    fn plain(path: &str) -> (String, String, bool) {
        (path.to_string(), "-".to_string(), false)
    }

    fn canonical(path: &str) -> (String, String, bool) {
        (path.to_string(), "-".to_string(), true)
    }

    #[test]
    fn the_path_is_the_positional_or_the_input_option() {
        assert_eq!(run(&["out.parquet"]), plain("out.parquet"));
        assert_eq!(run(&["--input", "a"]), plain("a"));
        assert_eq!(run(&["-i", "a"]), plain("a"));
        assert_eq!(run(&["--input=a"]), plain("a"));
        assert_eq!(run(&["-ia"]), plain("a"));
        assert_eq!(run(&["-i=a"]), plain("a"));
        assert_eq!(run(&["-ih"]), plain("h"));
        assert_eq!(run(&["--inp=a=b"]), plain("a=b"));
        assert_eq!(run(&["--input", "a", "b"]), plain("a"));
        assert_eq!(run(&["--input", "", "b"]), plain("b"));
        assert_eq!(run(&["--input", "no", "--input", "a"]), plain("a"));
        assert_eq!(run(&["--input=-x"]), plain("-x"));
        assert_eq!(run(&["--input", "-1"]), plain("-1"));
        assert_eq!(run(&["--input", "-"]), plain("-"));
    }

    #[test]
    fn options_take_every_spelling_argparse_takes() {
        assert_eq!(run(&["--canonical", "a"]), canonical("a"));
        assert_eq!(run(&["a", "--canonical", "--canonical"]), canonical("a"));
        assert_eq!(run(&["--c", "a"]), canonical("a"));
        assert_eq!(run(&["--canon", "a"]), canonical("a"));
        assert_eq!(run(&["--out", "-", "--in", "a", "--can"]), canonical("a"));
        let to = |output: &str| ("a".to_string(), output.to_string(), false);
        assert_eq!(run(&["a", "--output", "o"]), to("o"));
        assert_eq!(run(&["a", "-o", "o"]), to("o"));
        assert_eq!(run(&["a", "-oo"]), to("o"));
        assert_eq!(run(&["a", "-o-"]), to("-"));
        assert_eq!(run(&["a", "--o=o"]), to("o"));
        assert_eq!(run(&["a", "--output="]), to(""));
        assert_eq!(run(&["a", "-o", "first", "-o", "-"]), to("-"));
    }

    #[test]
    fn dashes_that_name_no_option_are_positional() {
        assert_eq!(run(&["--canonical", "--", "a"]), canonical("a"));
        assert_eq!(run(&["--", "-x"]), plain("-x"));
        assert_eq!(run(&["-"]), plain("-"));
        assert_eq!(run(&["-1"]), plain("-1"));
        assert_eq!(run(&["-1.5"]), plain("-1.5"));
        assert_eq!(run(&["-.5"]), plain("-.5"));
        assert_eq!(run(&["-a b"]), plain("-a b"));
        assert_eq!(run(&["--a b"]), plain("--a b"));
    }

    #[test]
    fn help_is_asked_where_it_stands() {
        assert_eq!(parsed(&["--help"]), Parsed::Help);
        assert_eq!(parsed(&["-h"]), Parsed::Help);
        assert_eq!(parsed(&["--he"]), Parsed::Help);
        assert_eq!(parsed(&["-hx"]), Parsed::Help);
        assert_eq!(parsed(&["--bogus", "--help"]), Parsed::Help);
        assert_eq!(parsed(&["a", "b", "--help"]), Parsed::Help);
        assert_eq!(
            usage(&["-i", "--help"]),
            "argument --input/-i: expected one argument"
        );
        assert_eq!(
            usage(&["-h=x"]),
            "argument -h/--help: ignored explicit argument 'x'"
        );
    }

    #[test]
    fn usage_errors_read_as_argparse_words_them() {
        assert_eq!(usage(&[]), "a Parquet path is required");
        assert_eq!(usage(&[""]), "a Parquet path is required");
        assert_eq!(usage(&["--canonical"]), "a Parquet path is required");
        assert_eq!(usage(&["", "a"]), "unrecognized arguments: a");
        assert_eq!(usage(&["a", "b"]), "unrecognized arguments: b");
        assert_eq!(
            usage(&["a", "--canonical", "b"]),
            "unrecognized arguments: b"
        );
        assert_eq!(usage(&["--", "--", "a"]), "unrecognized arguments: a");
        assert_eq!(usage(&["--bogus", "a"]), "unrecognized arguments: --bogus");
        assert_eq!(usage(&["-c", "a"]), "unrecognized arguments: -c");
        assert_eq!(
            usage(&["--canonicalx", "a"]),
            "unrecognized arguments: --canonicalx"
        );
        assert_eq!(
            usage(&["a", "--outputx", "-"]),
            "unrecognized arguments: --outputx -"
        );
        assert_eq!(usage(&["-i"]), "argument --input/-i: expected one argument");
        assert_eq!(
            usage(&["a", "-o"]),
            "argument --output/-o: expected one argument"
        );
        assert_eq!(
            usage(&["--input", "-x"]),
            "argument --input/-i: expected one argument"
        );
        assert_eq!(
            usage(&["--input", "--", "a"]),
            "argument --input/-i: expected one argument"
        );
        assert_eq!(
            usage(&["--canonical=x", "a"]),
            "argument --canonical: ignored explicit argument 'x'"
        );
        assert_eq!(
            usage(&["--=x", "a"]),
            "ambiguous option: --=x could match --help, --input, --output, --canonical, --types-of"
        );
    }

    #[test]
    fn the_types_option_names_a_parquet_output() {
        let named = |args: &[&str]| match parsed(args) {
            Parsed::Run(args) => args.types_of.map(|named| lossy(&named)),
            other => panic!("{args:?}: {other:?}"),
        };
        assert_eq!(named(&["--canonical", "a"]), None);
        assert_eq!(
            named(&["--canonical", "--types-of", "out.parquet", "a"]),
            Some("out.parquet".to_string())
        );
        assert_eq!(named(&["a", "--types-of=p"]), Some("p".to_string()));
        assert_eq!(named(&["a", "--t", "p"]), Some("p".to_string()));
        assert_eq!(
            named(&["a", "--types", "p", "--types-of", "q"]),
            Some("q".to_string())
        );
        assert_eq!(named(&["a", "--types-of", ""]), None);
        assert_eq!(run(&["--types-of", "p", "a"]), plain("a"));
        assert_eq!(
            usage(&["a", "--types-of"]),
            "argument --types-of: expected one argument"
        );
        assert_eq!(usage(&["a", "-t", "p"]), "unrecognized arguments: -t p");
    }
}
