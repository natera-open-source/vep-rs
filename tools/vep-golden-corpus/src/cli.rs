// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as CPython's `argparse` reads a parser of seven subcommands with
//! its defaults (`allow_abbrev`, options and positionals in any order): a long option is
//! named whole or by a prefix that fits it alone and takes its value from the next argument
//! or after `=`, a repeatable option keeps every value and any other its last one, and an
//! argument is acted on as it is read, so the first refusal or `-h` from the left ends the
//! reading. `generated` reads no abbreviated option, and the words after its `--command`
//! are set aside before the rest is read. The help texts are laid out for 80 columns, as
//! `argparse` lays them out where it finds no terminal width.

use crate::constants::{DEFAULT_INPUT, DEFAULT_MATRIX, PERL_MODES};
use crate::pytext::{py_int, py_repr_str};

/// The options of `select`. A repeatable option that is not given is an empty list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SelectArgs {
    pub release: String,
    pub assembly: String,
    pub species: Option<String>,
    pub dataset: Vec<String>,
    pub out: String,
    pub k: i64,
    pub seed: i64,
    pub max_records: i64,
    pub multi_allelic_per_dataset: i64,
    pub max_span: i64,
    pub focus_term: Vec<String>,
    pub focus_k: Option<i64>,
    pub include: Vec<String>,
    pub keep_records: bool,
    pub cache: Option<String>,
    pub flank: i64,
    pub covers: Vec<String>,
    pub family_columns: Option<String>,
    pub keep_samples: i64,
    pub genotype_exemplars: i64,
    pub require_genotype_class: Vec<String>,
    pub record_classes: Option<String>,
    pub record_class_exemplars: i64,
    pub matrix: Option<String>,
    pub drop_reference_statements: bool,
    pub flag: Vec<String>,
    pub fasta_name: Option<String>,
}

/// The options of `classify`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassifyArgs {
    pub corpus: String,
    pub vep_default: String,
    pub vep_rs_output: String,
    /// The reference's output of each further format (tab, VCF, JSON).
    pub vep_tab: Option<String>,
    pub vep_vcf: Option<String>,
    pub vep_json: Option<String>,
    /// vep-rs's output of each further format.
    pub vep_rs_tab: Option<String>,
    pub vep_rs_vcf: Option<String>,
    pub vep_rs_json: Option<String>,
    pub field: Vec<String>,
    pub flag_set: Option<String>,
    pub spliceai_table: Vec<String>,
    pub reference_fasta: Option<String>,
    pub formats_comparator: Option<String>,
    pub flag: Vec<String>,
    pub fasta_name: Option<String>,
}

/// The options of `check`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CheckArgs {
    pub corpus: String,
    pub dataset: Vec<String>,
    pub matrix: Option<String>,
    pub family_columns: Option<String>,
    pub require_genotype_class: Vec<String>,
}

/// The paths and the option of `identifiers`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdentifiersArgs {
    pub path: Vec<String>,
    pub check: bool,
}

/// The corpus and the options of `covers`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoversArgs {
    pub corpus: String,
    pub add: Vec<String>,
    pub list: bool,
    pub flag_set: Option<String>,
    pub matrix: String,
}

/// The corpus and the options of `flag-set`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlagSetArgs {
    pub corpus: String,
    pub name: Option<String>,
    pub flag: Vec<String>,
    pub input: String,
    pub stdin: bool,
    pub unordered: bool,
    pub perl_mode: String,
    pub format: Option<String>,
    pub expected_dir: Option<String>,
    pub covers: Vec<String>,
    pub matrix: String,
    pub list: bool,
    pub drop_reference_statements: bool,
}

/// The options of `generated`. `command` is `None` where no `--command` is given and holds
/// every word after the option where one is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneratedArgs {
    pub corpus: String,
    pub check: bool,
    pub command: Option<Vec<String>>,
}

/// One subcommand with its options.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Command {
    Select(Box<SelectArgs>),
    Classify(ClassifyArgs),
    Check(CheckArgs),
    Identifiers(IdentifiersArgs),
    Covers(CoversArgs),
    FlagSet(FlagSetArgs),
    Generated(GeneratedArgs),
}

/// What a command line asks for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Parsed {
    Run(Box<Command>),
    /// The help text for standard output; the run ends with exit code 0.
    Help(&'static str),
    /// The usage and the error line for standard error; the run ends with exit code 2.
    Refused(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Action {
    Help,
    Flag,
    Store,
    Append,
    /// Every argument after the option, options among them, up to a `--`.
    Remainder,
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    Text,
    Int,
    Choice(&'static [&'static str]),
}

#[derive(Debug)]
struct Opt {
    long: &'static str,
    short: Option<&'static str>,
    action: Action,
    kind: Kind,
    required: bool,
}

impl Opt {
    /// The option as an error names it: its strings joined by `/`.
    fn name(&self) -> String {
        match self.short {
            Some(short) => format!("{short}/{}", self.long),
            None => self.long.to_string(),
        }
    }

    fn takes_value(&self) -> bool {
        matches!(self.action, Action::Store | Action::Append)
    }
}

const fn option(long: &'static str, action: Action, kind: Kind, required: bool) -> Opt {
    Opt {
        long,
        short: None,
        action,
        kind,
        required,
    }
}

const fn store(long: &'static str) -> Opt {
    option(long, Action::Store, Kind::Text, false)
}

const fn required(long: &'static str) -> Opt {
    option(long, Action::Store, Kind::Text, true)
}

const fn int(long: &'static str) -> Opt {
    option(long, Action::Store, Kind::Int, false)
}

const fn append(long: &'static str) -> Opt {
    option(long, Action::Append, Kind::Text, false)
}

const fn flag(long: &'static str) -> Opt {
    option(long, Action::Flag, Kind::Text, false)
}

const HELP: Opt = Opt {
    long: "--help",
    short: Some("-h"),
    action: Action::Help,
    kind: Kind::Text,
    required: false,
};

#[derive(Debug)]
enum Positional {
    /// One argument.
    One(&'static str),
    /// One argument or more, up to the next option.
    OneOrMore(&'static str),
    /// A subcommand's name and every argument after it.
    Subcommand(&'static str),
}

impl Positional {
    fn name(&self) -> &'static str {
        match self {
            Positional::One(name) | Positional::OneOrMore(name) | Positional::Subcommand(name) => {
                name
            }
        }
    }
}

/// One parser: its options in the order they were declared, after its positional.
#[derive(Debug)]
struct Parser {
    prog: &'static str,
    help: &'static str,
    options: &'static [Opt],
    positional: Option<Positional>,
    /// Whether a prefix of a long option names it (`allow_abbrev`).
    abbreviations: bool,
}

impl Parser {
    /// The usage block of the help text, as an error prints it.
    fn usage(&self) -> &'static str {
        let end = self.help.find("\n\n").expect("a usage block") + 1;
        &self.help[..end]
    }

    /// Every option string with its option, in the order `argparse` holds them.
    fn option_strings(&'static self) -> impl Iterator<Item = (&'static str, &'static Opt)> {
        self.options.iter().flat_map(|option| {
            option
                .short
                .into_iter()
                .chain([option.long])
                .map(move |string| (string, option))
        })
    }

    fn exact(&'static self, string: &str) -> Option<&'static Opt> {
        self.option_strings()
            .find(|(held, _)| *held == string)
            .map(|(_, option)| option)
    }

    /// The two lines a refused command line ends on.
    fn refusal(&self, message: &str) -> String {
        format!("{}{}: error: {message}\n", self.usage(), self.prog)
    }
}

static TOP: Parser = Parser {
    prog: "vep-golden-corpus",
    help: include_str!("help/top.txt"),
    options: &[HELP],
    positional: Some(Positional::Subcommand("cmd")),
    abbreviations: true,
};

static SELECT: Parser = Parser {
    prog: "vep-golden-corpus select",
    help: include_str!("help/select.txt"),
    options: &[
        HELP,
        required("--release"),
        required("--assembly"),
        store("--species"),
        option("--dataset", Action::Append, Kind::Text, true),
        required("--out"),
        int("--k"),
        int("--seed"),
        int("--max-records"),
        int("--multi-allelic-per-dataset"),
        int("--max-span"),
        append("--focus-term"),
        int("--focus-k"),
        append("--include"),
        flag("--keep-records"),
        store("--cache"),
        int("--flank"),
        append("--covers"),
        store("--family-columns"),
        int("--keep-samples"),
        int("--genotype-exemplars"),
        append("--require-genotype-class"),
        store("--record-classes"),
        int("--record-class-exemplars"),
        store("--matrix"),
        flag("--drop-reference-statements"),
        append("--flag"),
        store("--fasta-name"),
    ],
    positional: None,
    abbreviations: true,
};

static CLASSIFY: Parser = Parser {
    prog: "vep-golden-corpus classify",
    help: include_str!("help/classify.txt"),
    options: &[
        HELP,
        required("--corpus"),
        required("--vep-default"),
        required("--vep-rs-output"),
        store("--vep-tab"),
        store("--vep-vcf"),
        store("--vep-json"),
        store("--vep-rs-tab"),
        store("--vep-rs-vcf"),
        store("--vep-rs-json"),
        append("--field"),
        store("--flag-set"),
        append("--spliceai-table"),
        store("--reference-fasta"),
        store("--formats-comparator"),
        append("--flag"),
        store("--fasta-name"),
    ],
    positional: None,
    abbreviations: true,
};

static CHECK: Parser = Parser {
    prog: "vep-golden-corpus check",
    help: include_str!("help/check.txt"),
    options: &[
        HELP,
        required("--corpus"),
        option("--dataset", Action::Append, Kind::Text, true),
        store("--matrix"),
        store("--family-columns"),
        append("--require-genotype-class"),
    ],
    positional: None,
    abbreviations: true,
};

static IDENTIFIERS: Parser = Parser {
    prog: "vep-golden-corpus identifiers",
    help: include_str!("help/identifiers.txt"),
    options: &[HELP, flag("--check")],
    positional: Some(Positional::OneOrMore("path")),
    abbreviations: true,
};

static COVERS: Parser = Parser {
    prog: "vep-golden-corpus covers",
    help: include_str!("help/covers.txt"),
    options: &[
        HELP,
        append("--add"),
        flag("--list"),
        store("--flag-set"),
        store("--matrix"),
    ],
    positional: Some(Positional::One("corpus")),
    abbreviations: true,
};

static FLAG_SET: Parser = Parser {
    prog: "vep-golden-corpus flag-set",
    help: include_str!("help/flag-set.txt"),
    options: &[
        HELP,
        store("--name"),
        append("--flag"),
        store("--input"),
        flag("--stdin"),
        flag("--unordered"),
        option(
            "--perl-mode",
            Action::Store,
            Kind::Choice(&PERL_MODES),
            false,
        ),
        store("--format"),
        store("--expected-dir"),
        append("--covers"),
        store("--matrix"),
        flag("--list"),
        flag("--drop-reference-statements"),
    ],
    positional: Some(Positional::One("corpus")),
    abbreviations: true,
};

/// The option of `generated` whose words are the command.
const COMMAND_OPTION: &str = "--command";

static GENERATED: Parser = Parser {
    prog: "vep-golden-corpus generated",
    help: include_str!("help/generated.txt"),
    options: &[
        HELP,
        required("--corpus"),
        flag("--check"),
        option(COMMAND_OPTION, Action::Remainder, Kind::Text, false),
    ],
    positional: None,
    abbreviations: false,
};

/// The name of the subcommand whose `--command` takes every word after it.
const GENERATED_NAME: &str = "generated";

/// The subcommands in the order an error lists them.
const SUBCOMMANDS: [(&str, &Parser); 7] = [
    ("select", &SELECT),
    ("classify", &CLASSIFY),
    ("check", &CHECK),
    ("identifiers", &IDENTIFIERS),
    ("covers", &COVERS),
    ("flag-set", &FLAG_SET),
    (GENERATED_NAME, &GENERATED),
];

/// How one argument that starts like an option may be read.
#[derive(Clone, Debug)]
struct Candidate<'a> {
    /// The option, `None` for an argument shaped like an option that names none.
    option: Option<&'static Opt>,
    /// The option string it was matched under.
    string: String,
    /// What stands between the option string and its explicit argument: `=`, or nothing for
    /// the characters behind a one-letter option.
    separator: Option<&'static str>,
    explicit: Option<&'a str>,
}

/// What the reading hands on, in the order of the command line.
enum Taken<'a> {
    /// An option, with its value where it takes one.
    Option(&'static Opt, Option<&'a str>),
    /// The arguments the option that takes the remainder took.
    Words(Vec<&'a str>),
    /// The arguments a positional takes.
    Positional(Vec<&'a str>),
}

/// What ends a reading before its last argument.
enum Stop {
    Help(&'static str),
    /// The whole text for standard error.
    Refused(String),
    /// An error of the parser in hand, without its usage.
    Error(String),
}

/// Whether an argument that starts with `-` reads as an option, and as which.
fn parse_optional<'a>(parser: &'static Parser, arg: &'a str) -> Option<Vec<Candidate<'a>>> {
    if !arg.starts_with('-') {
        return None;
    }
    let known = |option, string: &str, separator, explicit| {
        Some(vec![Candidate {
            option,
            string: string.to_string(),
            separator,
            explicit,
        }])
    };
    if let Some(option) = parser.exact(arg) {
        return known(Some(option), arg, None, None);
    }
    if arg == "-" {
        return None;
    }
    if let Some((head, explicit)) = arg.split_once('=') {
        if let Some(option) = parser.exact(head) {
            return known(Some(option), head, Some("="), Some(explicit));
        }
    }
    let candidates = option_tuples(parser, arg);
    if !candidates.is_empty() {
        return Some(candidates);
    }
    // An argument that starts like a negative number, or holds a space, is a value.
    let digits = arg[1..].strip_prefix('.').unwrap_or(&arg[1..]);
    if digits.starts_with(|c: char| c.is_ascii_digit()) || arg.contains(' ') {
        return None;
    }
    known(None, arg, None, None)
}

/// The options an argument abbreviates: every long option its text before `=` starts, in a
/// parser that reads abbreviations, and for a one-dash argument the one-letter option of
/// its first two characters with the rest as that option's explicit argument.
fn option_tuples<'a>(parser: &'static Parser, arg: &'a str) -> Vec<Candidate<'a>> {
    let (prefix, separator, explicit) = match arg.split_once('=') {
        Some((head, tail)) => (head, Some("="), Some(tail)),
        None => (arg, None, None),
    };
    let second = arg[1..].chars().next().expect("two characters or more");
    let short_end = 1 + second.len_utf8();
    let mut candidates = Vec::new();
    for (string, option) in parser.option_strings() {
        if second != '-' && string == &arg[..short_end] {
            candidates.push(Candidate {
                option: Some(option),
                string: string.to_string(),
                separator: Some(""),
                explicit: Some(&arg[short_end..]),
            });
        } else if parser.abbreviations && string.starts_with(prefix) {
            candidates.push(Candidate {
                option: Some(option),
                string: string.to_string(),
                separator,
                explicit,
            });
        }
    }
    candidates
}

/// One reading of an argument list by one parser.
struct Reading<'a, 'f> {
    parser: &'static Parser,
    args: &'a [&'a str],
    /// One byte an argument: `A` for a value, `O` for an option, `-` for the `--` marker.
    pattern: Vec<u8>,
    candidates: Vec<Option<Vec<Candidate<'a>>>>,
    extras: Vec<String>,
    seen: Vec<&'static str>,
    positional_left: bool,
    take: &'f mut dyn FnMut(Taken<'a>) -> Result<(), Stop>,
}

impl<'a> Reading<'a, '_> {
    /// Reads the option at `start` and returns the index after it and its value.
    fn consume_optional(&mut self, start: usize) -> Result<usize, Stop> {
        let candidates = self.candidates[start].clone().expect("an option");
        if candidates.len() > 1 {
            let matches: Vec<&str> = candidates.iter().map(|c| c.string.as_str()).collect();
            return Err(Stop::Error(format!(
                "ambiguous option: {} could match {}",
                self.args[start],
                matches.join(", ")
            )));
        }
        let Candidate {
            mut option,
            mut string,
            mut separator,
            mut explicit,
        } = candidates.into_iter().next().expect("one candidate");
        let mut taken: Vec<(&'static Opt, Option<&'a str>)> = Vec::new();
        let stop;
        loop {
            let Some(found) = option else {
                self.extras.push(self.args[start].to_string());
                return Ok(start + 1);
            };
            if found.action == Action::Remainder {
                // The option's own argument after `=` is its one word; without one it takes
                // every argument up to a `--`, which no option reads past.
                let (words, after) = match explicit {
                    Some(word) => (vec![word], start + 1),
                    None => {
                        let rest = &self.pattern[start + 1..];
                        let count = rest.iter().take_while(|byte| **byte != b'-').count();
                        let end = start + 1 + count;
                        (self.args[start + 1..end].to_vec(), end)
                    }
                };
                self.seen.push(found.long);
                (self.take)(Taken::Words(words))?;
                return Ok(after);
            }
            let ignored = |explicit: &str| {
                Stop::Error(format!(
                    "argument {}: ignored explicit argument {}",
                    found.name(),
                    py_repr_str(explicit)
                ))
            };
            let Some(argument) = explicit else {
                let values = usize::from(found.takes_value());
                if values == 1 && self.pattern.get(start + 1) != Some(&b'A') {
                    return Err(Stop::Error(format!(
                        "argument {}: expected one argument",
                        found.name()
                    )));
                }
                taken.push((found, (values == 1).then(|| self.args[start + 1])));
                stop = start + 1 + values;
                break;
            };
            if found.takes_value() {
                taken.push((found, Some(argument)));
                stop = start + 1;
                break;
            }
            // The characters behind a one-letter option that takes no value read as more
            // one-letter options.
            let one_letter = !string[1..].starts_with('-');
            if !one_letter || argument.is_empty() {
                return Err(ignored(argument));
            }
            if separator == Some("=") || argument.starts_with('-') {
                return Err(ignored(argument));
            }
            taken.push((found, None));
            let letter = argument.chars().next().expect("a character");
            let rest = &argument[letter.len_utf8()..];
            string = format!("-{letter}");
            match self.parser.exact(&string) {
                Some(next) => {
                    option = Some(next);
                    (separator, explicit) = match rest.strip_prefix('=') {
                        _ if rest.is_empty() => (None, None),
                        Some(value) => (Some("="), Some(value)),
                        None => (Some(""), Some(rest)),
                    };
                }
                None => {
                    self.extras.push(format!("-{argument}"));
                    stop = start + 1;
                    break;
                }
            }
        }
        for (found, value) in taken {
            self.seen.push(found.long);
            (self.take)(Taken::Option(found, value))?;
        }
        Ok(stop)
    }

    /// Hands the positional the arguments it takes from `start` and returns the index after
    /// them; `start` where it takes none.
    fn consume_positionals(&mut self, start: usize) -> Result<usize, Stop> {
        let Some(positional) = self
            .parser
            .positional
            .as_ref()
            .filter(|_| self.positional_left)
        else {
            return Ok(start);
        };
        let pattern = &self.pattern[start..];
        let markers = pattern.iter().take_while(|byte| **byte == b'-').count();
        if pattern.get(markers) != Some(&b'A') {
            return Ok(start);
        }
        let count = match positional {
            Positional::Subcommand(_) => pattern.len(),
            Positional::One(_) => {
                let trailing = pattern[markers + 1..]
                    .iter()
                    .take_while(|byte| **byte == b'-')
                    .count();
                markers + 1 + trailing
            }
            Positional::OneOrMore(_) => {
                let further = pattern[markers + 1..]
                    .iter()
                    .take_while(|byte| **byte != b'O')
                    .count();
                markers + 1 + further
            }
        };
        let mut values: Vec<&'a str> = self.args[start..start + count].to_vec();
        let marker = match positional {
            Positional::Subcommand(_) => (markers > 0).then_some(0),
            Positional::One(_) | Positional::OneOrMore(_) => {
                pattern[..count].iter().position(|byte| *byte == b'-')
            }
        };
        if let Some(place) = marker {
            values.remove(place);
        }
        self.positional_left = false;
        self.seen.push(positional.name());
        (self.take)(Taken::Positional(values))?;
        Ok(start + count)
    }
}

/// Reads `args` with `parser`, handing each option and positional to `take` as it is read,
/// and returns the arguments nothing took.
fn parse_known<'a>(
    parser: &'static Parser,
    args: &'a [&'a str],
    take: &mut dyn FnMut(Taken<'a>) -> Result<(), Stop>,
) -> Result<Vec<String>, Stop> {
    let mut pattern = Vec::with_capacity(args.len());
    let mut candidates = Vec::with_capacity(args.len());
    let mut values_only = false;
    for arg in args {
        let read = match *arg {
            _ if values_only => None,
            "--" => {
                values_only = true;
                pattern.push(b'-');
                candidates.push(None);
                continue;
            }
            _ => parse_optional(parser, arg),
        };
        pattern.push(if read.is_some() { b'O' } else { b'A' });
        candidates.push(read);
    }
    let last_option = candidates.iter().rposition(Option::is_some);
    let mut reading = Reading {
        parser,
        args,
        pattern,
        candidates,
        extras: Vec::new(),
        seen: Vec::new(),
        positional_left: parser.positional.is_some(),
        take,
    };
    let mut start = 0;
    while let Some(last_option) = last_option.filter(|last| start <= *last) {
        let next_option = (start..=last_option)
            .find(|index| reading.candidates[*index].is_some())
            .expect("an option at or before the last one");
        if start != next_option {
            let end = reading.consume_positionals(start)?;
            if end > start {
                start = end;
                continue;
            }
        }
        if reading.candidates[start].is_none() {
            let unread = args[start..next_option].iter().map(|arg| arg.to_string());
            reading.extras.extend(unread);
            start = next_option;
        }
        start = reading.consume_optional(start)?;
    }
    let stop = reading.consume_positionals(start)?;
    reading
        .extras
        .extend(args[stop..].iter().map(|arg| arg.to_string()));
    let positional = parser.positional.iter().map(Positional::name);
    let options = parser
        .options
        .iter()
        .filter(|option| option.required)
        .map(|option| option.long);
    let missing: Vec<&str> = positional
        .chain(options)
        .filter(|name| !reading.seen.contains(name))
        .collect();
    if !missing.is_empty() {
        return Err(Stop::Error(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    Ok(reading.extras)
}

/// What a subcommand's reading holds: every option's values in the order given.
#[derive(Default)]
struct Matches {
    values: Vec<(&'static str, String)>,
    flags: Vec<&'static str>,
    positional: Vec<String>,
    /// The words of the option that takes the remainder, as it was given last.
    words: Option<Vec<String>>,
}

impl Matches {
    fn take(&mut self, parser: &'static Parser, taken: Taken<'_>) -> Result<(), Stop> {
        let (option, value) = match taken {
            Taken::Positional(values) => {
                self.positional = values.iter().map(|value| value.to_string()).collect();
                return Ok(());
            }
            Taken::Words(words) => {
                self.words = Some(words.iter().map(|word| word.to_string()).collect());
                return Ok(());
            }
            Taken::Option(option, value) => (option, value),
        };
        match option.action {
            Action::Help => return Err(Stop::Help(parser.help)),
            Action::Flag => {
                self.flags.push(option.long);
                return Ok(());
            }
            Action::Store | Action::Append | Action::Remainder => {}
        }
        let value = value.expect("a value for an option that takes one");
        let refused = match option.kind {
            Kind::Text => None,
            Kind::Int => py_int(value)
                .is_none()
                .then(|| format!("invalid int value: {}", py_repr_str(value))),
            Kind::Choice(choices) => (!choices.contains(&value)).then(|| {
                let listed: Vec<String> = choices.iter().map(|c| py_repr_str(c)).collect();
                format!(
                    "invalid choice: {} (choose from {})",
                    py_repr_str(value),
                    listed.join(", ")
                )
            }),
        };
        if let Some(message) = refused {
            return Err(Stop::Error(format!(
                "argument {}: {message}",
                option.name()
            )));
        }
        self.values.push((option.long, value.to_string()));
        Ok(())
    }

    fn last(&self, name: &str) -> Option<String> {
        let found = self.values.iter().rev().find(|(held, _)| *held == name);
        found.map(|(_, value)| value.clone())
    }

    fn all(&self, name: &str) -> Vec<String> {
        let found = self.values.iter().filter(|(held, _)| *held == name);
        found.map(|(_, value)| value.clone()).collect()
    }

    fn flag(&self, name: &str) -> bool {
        self.flags.contains(&name)
    }

    fn int(&self, name: &str) -> Option<i64> {
        self.last(name)
            .map(|value| py_int(&value).expect("a count held as it was read"))
    }

    fn required(&self, name: &str) -> String {
        self.last(name).expect("a required option")
    }

    fn corpus(&self) -> String {
        self.positional.first().cloned().expect("the corpus")
    }

    fn command(&self, name: &str) -> Command {
        match name {
            "select" => Command::Select(Box::new(SelectArgs {
                release: self.required("--release"),
                assembly: self.required("--assembly"),
                species: self.last("--species"),
                dataset: self.all("--dataset"),
                out: self.required("--out"),
                k: self.int("--k").unwrap_or(3),
                seed: self.int("--seed").unwrap_or(1),
                max_records: self.int("--max-records").unwrap_or(1500),
                multi_allelic_per_dataset: self.int("--multi-allelic-per-dataset").unwrap_or(4),
                max_span: self.int("--max-span").unwrap_or(1_000_000),
                focus_term: self.all("--focus-term"),
                focus_k: self.int("--focus-k"),
                include: self.all("--include"),
                keep_records: self.flag("--keep-records"),
                cache: self.last("--cache"),
                flank: self.int("--flank").unwrap_or(5000),
                covers: self.all("--covers"),
                family_columns: self.last("--family-columns"),
                keep_samples: self.int("--keep-samples").unwrap_or(0),
                genotype_exemplars: self.int("--genotype-exemplars").unwrap_or(8),
                require_genotype_class: self.all("--require-genotype-class"),
                record_classes: self.last("--record-classes"),
                record_class_exemplars: self.int("--record-class-exemplars").unwrap_or(2),
                matrix: self.last("--matrix"),
                drop_reference_statements: self.flag("--drop-reference-statements"),
                flag: self.all("--flag"),
                fasta_name: self.last("--fasta-name"),
            })),
            "classify" => Command::Classify(ClassifyArgs {
                corpus: self.required("--corpus"),
                vep_default: self.required("--vep-default"),
                vep_rs_output: self.required("--vep-rs-output"),
                vep_tab: self.last("--vep-tab"),
                vep_vcf: self.last("--vep-vcf"),
                vep_json: self.last("--vep-json"),
                vep_rs_tab: self.last("--vep-rs-tab"),
                vep_rs_vcf: self.last("--vep-rs-vcf"),
                vep_rs_json: self.last("--vep-rs-json"),
                field: self.all("--field"),
                flag_set: self.last("--flag-set"),
                spliceai_table: self.all("--spliceai-table"),
                reference_fasta: self.last("--reference-fasta"),
                formats_comparator: self.last("--formats-comparator"),
                flag: self.all("--flag"),
                fasta_name: self.last("--fasta-name"),
            }),
            "check" => Command::Check(CheckArgs {
                corpus: self.required("--corpus"),
                dataset: self.all("--dataset"),
                matrix: self.last("--matrix"),
                family_columns: self.last("--family-columns"),
                require_genotype_class: self.all("--require-genotype-class"),
            }),
            "identifiers" => Command::Identifiers(IdentifiersArgs {
                path: self.positional.clone(),
                check: self.flag("--check"),
            }),
            "covers" => Command::Covers(CoversArgs {
                corpus: self.corpus(),
                add: self.all("--add"),
                list: self.flag("--list"),
                flag_set: self.last("--flag-set"),
                matrix: self
                    .last("--matrix")
                    .unwrap_or_else(|| DEFAULT_MATRIX.to_string()),
            }),
            GENERATED_NAME => Command::Generated(GeneratedArgs {
                corpus: self.required("--corpus"),
                check: self.flag("--check"),
                command: self.words.clone(),
            }),
            _ => Command::FlagSet(FlagSetArgs {
                corpus: self.corpus(),
                name: self.last("--name"),
                flag: self.all("--flag"),
                input: self
                    .last("--input")
                    .unwrap_or_else(|| DEFAULT_INPUT.to_string()),
                stdin: self.flag("--stdin"),
                unordered: self.flag("--unordered"),
                perl_mode: self
                    .last("--perl-mode")
                    .unwrap_or_else(|| "offline".to_string()),
                format: self.last("--format"),
                expected_dir: self.last("--expected-dir"),
                covers: self.all("--covers"),
                matrix: self
                    .last("--matrix")
                    .unwrap_or_else(|| DEFAULT_MATRIX.to_string()),
                list: self.flag("--list"),
                drop_reference_statements: self.flag("--drop-reference-statements"),
            }),
        }
    }
}

/// Reads the arguments after the program's name. Where the first is `generated`, every
/// argument after the first `--command` is a word of the command, a `--` and an option of
/// the subcommand among them, and the arguments before it are read alone.
pub fn parse(args: &[String]) -> Parsed {
    let cut = args
        .iter()
        .position(|arg| arg == COMMAND_OPTION)
        .filter(|_| args.first().is_some_and(|first| first == GENERATED_NAME));
    let Some(cut) = cut else {
        return parse_whole(args);
    };
    match parse_whole(&args[..cut]) {
        Parsed::Run(mut command) => {
            if let Command::Generated(generated) = command.as_mut() {
                generated.command = Some(args[cut + 1..].to_vec());
            }
            Parsed::Run(command)
        }
        ended => ended,
    }
}

/// Reads `args` as one command line.
fn parse_whole(args: &[String]) -> Parsed {
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let mut command = None;
    let mut unread_by_subcommand = Vec::new();
    let read = parse_known(&TOP, &args, &mut |taken| {
        let Taken::Positional(values) = taken else {
            return Err(Stop::Help(TOP.help));
        };
        let name = values[0];
        let Some((_, parser)) = SUBCOMMANDS.iter().find(|(held, _)| *held == name) else {
            let listed: Vec<String> = SUBCOMMANDS.iter().map(|(n, _)| py_repr_str(n)).collect();
            return Err(Stop::Error(format!(
                "argument cmd: invalid choice: {} (choose from {})",
                py_repr_str(name),
                listed.join(", ")
            )));
        };
        let mut matches = Matches::default();
        let unread = parse_known(parser, &values[1..], &mut |taken| {
            matches.take(parser, taken)
        })
        .map_err(|stop| match stop {
            Stop::Error(message) => Stop::Refused(parser.refusal(&message)),
            other => other,
        })?;
        unread_by_subcommand = unread;
        command = Some(matches.command(name));
        Ok(())
    });
    match read {
        Err(Stop::Help(text)) => Parsed::Help(text),
        Err(Stop::Refused(text)) => Parsed::Refused(text),
        Err(Stop::Error(message)) => Parsed::Refused(TOP.refusal(&message)),
        Ok(mut unread) => {
            unread.extend(unread_by_subcommand);
            if !unread.is_empty() {
                let message = format!("unrecognized arguments: {}", unread.join(" "));
                return Parsed::Refused(TOP.refusal(&message));
            }
            Parsed::Run(Box::new(
                command.expect("a subcommand, which the reading requires"),
            ))
        }
    }
}

#[cfg(test)]
#[path = "cli_parse_tests.rs"]
mod parse_tests;
