// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line of `score`, read as CPython's `argparse` reads one: an option by its
//! name, by `--name=value` or by a prefix that only one option has; a value that opens with
//! `-` taken only where it holds a space or reads as a negative number; the last value of
//! an option given twice, or every value of one that repeats; and `argparse`'s message for
//! a line it refuses.

use super::paths::path_text;
use super::text::{integer_text, is_decimal_digit, repr};

pub const USAGE: &str = "usage: vep-plugin-pairs score --reference OUTPUT --engine ID=PATH [--engine ID=PATH ...] --report JSON [--plugins A,B] [--custom-source NAME=PREFIX ...] [options]";
pub const HELP: &str = "\
Scores the plugin tuples of each engine's output against one Ensembl VEP default-format
output: per source, engine and projection the tuples of each side, their intersection, the
tuples a documented class sets aside and the residual rows.

options:
  -h, --help                      print this text
  --reference OUTPUT              Ensembl VEP default-format output, plain or gzip
  --engine ID=PATH                a scored engine's default-format output or tuple file; repeat per engine
  --engine-binary ID=SHA256       the sha256 of the binary that wrote an engine's output, carried into the report
  --engine-commit ID=COMMIT       the source commit of that binary, carried into the report
  --plugins A,B                   the plugins the scored run loaded, by their registry names or aliases
  --run-command COMMAND           the scored run's command line, whose flags the registry reads and whose
                                  assembly dbnsfp_mitochondrial_record_named_m reads
  --custom-source NAME=PREFIX     group the fields starting with PREFIX as the source NAME; repeat per source
  --work-dir DIR                  directory the scratch is made in (default: the system's temporary directory)
  --report JSON                   the report's path
  --residuals DIR                 directory for the residual files, <source>.<engine>.<projection>.tsv
  --class-table TSV               the documented classes: per line a class, its plugins and their field names
                                  (default: the table built into the program, plugin_class_table.tsv)
  --classes A,B                   the declared classes to apply (default: every one whose fourth column
                                  does not read 'not met'; an empty value: none)
  --class-population CLASS=N      the keys a class can reach in this input; a count above it is marked
  --input FILE                    the VCF the scored run annotated, read for the records its rows name
  --gwas-catalog FILE             the associations file the run's GWAS plugin read
  --revel-table FILE              the table the run's REVEL plugin read, its tabix index beside it
  --revel-no-match                the run set the REVEL plugin's no_match
  --utrannotator-sweep PATH       an Ensembl VEP output of the input under one hash seed, or a directory
                                  of them; repeat per run
  --json-cache DIR                the candidate JSON cache of the scored run
  --spliceai-table FILE           a score file the run's SpliceAI plugin read; repeat per file
  --pairing-threads N             threads of the run (default: every core)
  --pairing-tool BINARY           accepted for the command lines that name one: a file, not read
  --sort-memory SIZE              accepted and not read
  --sort-threads N                accepted and not read";

/// A whole number as `int` reads an option's value and `str` writes it back.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Integer(pub String);

/// The options of one command line, each as `argparse` leaves it: a path as
/// `str(pathlib.Path(value))` writes it, a repeated option's values in their order.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Args {
    pub reference: String,
    pub engine: Vec<String>,
    pub engine_binary: Vec<String>,
    pub engine_commit: Vec<String>,
    pub plugins: String,
    pub run_command: Option<String>,
    pub custom_source: Vec<String>,
    pub work_dir: Option<String>,
    pub report: String,
    pub residuals: Option<String>,
    pub class_table: Option<String>,
    pub gwas_catalog: Option<String>,
    pub revel_table: Option<String>,
    pub revel_no_match: bool,
    pub utrannotator_sweep: Vec<String>,
    pub json_cache: Option<String>,
    pub spliceai_table: Vec<String>,
    pub input: Option<String>,
    pub classes: Option<String>,
    pub class_population: Vec<String>,
    pub pairing_tool: Option<String>,
    pub pairing_threads: Option<Integer>,
    pub sort_memory: String,
    pub sort_threads: Option<Integer>,
}

/// What a command line asks for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Parsed {
    Help,
    Run(Box<Args>),
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Action {
    Help,
    Reference,
    Engine,
    EngineBinary,
    EngineCommit,
    Plugins,
    RunCommand,
    CustomSource,
    WorkDir,
    Report,
    Residuals,
    ClassTable,
    GwasCatalog,
    RevelTable,
    RevelNoMatch,
    UtrannotatorSweep,
    JsonCache,
    SpliceaiTable,
    Input,
    Classes,
    ClassPopulation,
    PairingTool,
    PairingThreads,
    SortMemory,
    SortThreads,
}

impl Action {
    fn takes_value(self) -> bool {
        !matches!(self, Action::Help | Action::RevelNoMatch)
    }

    /// The action's name in a message: its option strings joined by `/`.
    fn name(self) -> String {
        OPTIONS
            .iter()
            .filter(|(_, action)| *action == self)
            .map(|(option, _)| *option)
            .collect::<Vec<_>>()
            .join("/")
    }
}

/// Every option string in the order `argparse` holds them, which is the order it lists the
/// options a prefix could mean.
const OPTIONS: &[(&str, Action)] = &[
    ("-h", Action::Help),
    ("--help", Action::Help),
    ("--reference", Action::Reference),
    ("--engine", Action::Engine),
    ("--engine-binary", Action::EngineBinary),
    ("--engine-commit", Action::EngineCommit),
    ("--plugins", Action::Plugins),
    ("--run-command", Action::RunCommand),
    ("--custom-source", Action::CustomSource),
    ("--work-dir", Action::WorkDir),
    ("--report", Action::Report),
    ("--residuals", Action::Residuals),
    ("--class-table", Action::ClassTable),
    ("--gwas-catalog", Action::GwasCatalog),
    ("--revel-table", Action::RevelTable),
    ("--revel-no-match", Action::RevelNoMatch),
    ("--utrannotator-sweep", Action::UtrannotatorSweep),
    ("--json-cache", Action::JsonCache),
    ("--spliceai-table", Action::SpliceaiTable),
    ("--input", Action::Input),
    ("--classes", Action::Classes),
    ("--class-population", Action::ClassPopulation),
    ("--pairing-tool", Action::PairingTool),
    ("--pairing-threads", Action::PairingThreads),
    ("--sort-memory", Action::SortMemory),
    ("--sort-threads", Action::SortThreads),
];
const REQUIRED: [Action; 3] = [Action::Reference, Action::Engine, Action::Report];

fn action_of(option: &str) -> Option<Action> {
    OPTIONS
        .iter()
        .find(|(known, _)| *known == option)
        .map(|(_, action)| *action)
}

/// One reading of an argument as an option: the action, `None` for an option the program
/// lacks; the option string; and the value written in the same argument, with whether an
/// `=` stood before it.
#[derive(Clone, Debug)]
struct Reading {
    action: Option<Action>,
    option: String,
    separated: bool,
    explicit: Option<String>,
}

/// `str.partition("=")` with whether the separator was found.
fn at_equals(argument: &str) -> (&str, Option<&str>) {
    match argument.split_once('=') {
        Some((before, after)) => (before, Some(after)),
        None => (argument, None),
    }
}

/// The options `argument` could mean by a prefix, in the options' order.
fn by_prefix(argument: &str) -> Vec<Reading> {
    let mut characters = argument.chars();
    let (first, second) = (characters.next(), characters.next());
    let (prefix, explicit) = at_equals(argument);
    let mut found = Vec::new();
    if first == Some('-') && second == Some('-') {
        for (option, action) in OPTIONS {
            if option.starts_with(prefix) {
                found.push(Reading {
                    action: Some(*action),
                    option: option.to_string(),
                    separated: explicit.is_some(),
                    explicit: explicit.map(str::to_string),
                });
            }
        }
    } else {
        // A one-letter option takes the rest of its argument as its value.
        let letter_end = argument
            .char_indices()
            .nth(2)
            .map_or(argument.len(), |(at, _)| at);
        let (letter, rest) = argument.split_at(letter_end);
        for (option, action) in OPTIONS {
            if *option == letter {
                found.push(Reading {
                    action: Some(*action),
                    option: option.to_string(),
                    separated: false,
                    explicit: Some(rest.to_string()),
                });
            } else if option.starts_with(prefix) {
                found.push(Reading {
                    action: Some(*action),
                    option: option.to_string(),
                    separated: explicit.is_some(),
                    explicit: explicit.map(str::to_string),
                });
            }
        }
    }
    found
}

/// Whether the argument reads as a negative number: `-`, an optional `.`, a decimal digit.
fn negative_number(argument: &str) -> bool {
    let mut rest = argument.chars().skip(1);
    match rest.next() {
        Some('.') => rest.next().is_some_and(is_decimal_digit),
        Some(character) => is_decimal_digit(character),
        None => false,
    }
}

/// The readings of an argument as an option; `None` for a value.
fn as_option(argument: &str) -> Option<Vec<Reading>> {
    if !argument.starts_with('-') {
        return None;
    }
    if let Some(action) = action_of(argument) {
        return Some(vec![Reading {
            action: Some(action),
            option: argument.to_string(),
            separated: false,
            explicit: None,
        }]);
    }
    if argument.chars().count() == 1 {
        return None;
    }
    if let (option, Some(explicit)) = at_equals(argument) {
        if let Some(action) = action_of(option) {
            return Some(vec![Reading {
                action: Some(action),
                option: option.to_string(),
                separated: true,
                explicit: Some(explicit.to_string()),
            }]);
        }
    }
    let found = by_prefix(argument);
    if !found.is_empty() {
        return Some(found);
    }
    if negative_number(argument) || argument.contains(' ') {
        return None;
    }
    Some(vec![Reading {
        action: None,
        option: argument.to_string(),
        separated: false,
        explicit: None,
    }])
}

fn refusal(action: Action, message: &str) -> String {
    format!("argument {}: {message}", action.name())
}

fn whole_number(action: Action, value: &str) -> Result<Integer, String> {
    integer_text(value)
        .map(Integer)
        .ok_or_else(|| refusal(action, &format!("invalid int value: {}", repr(value))))
}

/// Takes one action with its value. `Ok(true)` ends the reading: the help was asked for.
fn take(args: &mut Args, action: Action, value: Option<&str>) -> Result<bool, String> {
    let text = || value.unwrap_or_default().to_string();
    let path = || path_text(value.unwrap_or_default());
    match action {
        Action::Help => return Ok(true),
        Action::Reference => args.reference = path(),
        Action::Engine => args.engine.push(text()),
        Action::EngineBinary => args.engine_binary.push(text()),
        Action::EngineCommit => args.engine_commit.push(text()),
        Action::Plugins => args.plugins = text(),
        Action::RunCommand => args.run_command = Some(text()),
        Action::CustomSource => args.custom_source.push(text()),
        Action::WorkDir => args.work_dir = Some(path()),
        Action::Report => args.report = path(),
        Action::Residuals => args.residuals = Some(path()),
        Action::ClassTable => args.class_table = Some(path()),
        Action::GwasCatalog => args.gwas_catalog = Some(path()),
        Action::RevelTable => args.revel_table = Some(path()),
        Action::RevelNoMatch => args.revel_no_match = true,
        Action::UtrannotatorSweep => args.utrannotator_sweep.push(path()),
        Action::JsonCache => args.json_cache = Some(path()),
        Action::SpliceaiTable => args.spliceai_table.push(path()),
        Action::Input => args.input = Some(path()),
        Action::Classes => args.classes = Some(text()),
        Action::ClassPopulation => args.class_population.push(text()),
        Action::PairingTool => args.pairing_tool = Some(path()),
        Action::PairingThreads => {
            args.pairing_threads = Some(whole_number(action, value.unwrap_or_default())?)
        }
        Action::SortMemory => args.sort_memory = text(),
        Action::SortThreads => {
            args.sort_threads = Some(whole_number(action, value.unwrap_or_default())?)
        }
    }
    Ok(false)
}

/// Reads the arguments after the program's word. The refusal is `argparse`'s message, the
/// text after `error: `.
pub fn parse(arguments: &[String]) -> Result<Parsed, String> {
    // Every argument is an option or a value before any is taken; after `--` all are values.
    let end = arguments
        .iter()
        .position(|argument| argument == "--")
        .unwrap_or(arguments.len());
    let readings: Vec<Option<Vec<Reading>>> = arguments
        .iter()
        .enumerate()
        .map(|(index, argument)| {
            if index < end {
                as_option(argument)
            } else {
                None
            }
        })
        .collect();
    let is_value = |index: usize| index < end && readings[index].is_none();

    let mut args = Args {
        sort_memory: "2G".to_string(),
        ..Args::default()
    };
    let mut seen: Vec<Action> = Vec::new();
    let mut unread: Vec<String> = Vec::new();
    let mut index = 0;
    while index < arguments.len() {
        let Some(found) = &readings[index] else {
            unread.push(arguments[index].clone());
            index += 1;
            continue;
        };
        if found.len() > 1 {
            let options: Vec<&str> = found
                .iter()
                .map(|reading| reading.option.as_str())
                .collect();
            return Err(format!(
                "ambiguous option: {} could match {}",
                arguments[index],
                options.join(", ")
            ));
        }
        let Reading {
            action,
            mut option,
            mut separated,
            mut explicit,
        } = found[0].clone();
        let Some(mut action) = action else {
            unread.push(arguments[index].clone());
            index += 1;
            continue;
        };
        // The actions one argument names, each with its value: a one-letter option that
        // takes none can stand before another in the same argument.
        let mut taken: Vec<(Action, Option<String>)> = Vec::new();
        let mut next = index + 1;
        loop {
            let Some(written) = explicit.clone() else {
                let value = if action.takes_value() {
                    if !is_value(next) {
                        return Err(refusal(action, "expected one argument"));
                    }
                    next += 1;
                    Some(arguments[next - 1].clone())
                } else {
                    None
                };
                taken.push((action, value));
                break;
            };
            if action.takes_value() {
                taken.push((action, Some(written)));
                break;
            }
            let one_letter = !option[1..].starts_with('-');
            if !one_letter || written.is_empty() || separated || written.starts_with('-') {
                return Err(refusal(
                    action,
                    &format!("ignored explicit argument {}", repr(&written)),
                ));
            }
            taken.push((action, None));
            let mut characters = written.chars();
            let letter = characters.next().unwrap_or_default();
            let rest = characters.as_str();
            option = format!("-{letter}");
            match action_of(&option) {
                Some(following) => {
                    action = following;
                    (separated, explicit) = match rest.strip_prefix('=') {
                        _ if rest.is_empty() => (false, None),
                        Some(value) => (true, Some(value.to_string())),
                        None => (false, Some(rest.to_string())),
                    };
                }
                None => {
                    unread.push(format!("-{written}"));
                    break;
                }
            }
        }
        for (action, value) in taken {
            if !seen.contains(&action) {
                seen.push(action);
            }
            if take(&mut args, action, value.as_deref())? {
                return Ok(Parsed::Help);
            }
        }
        index = next;
    }
    let missing: Vec<String> = REQUIRED
        .iter()
        .filter(|action| !seen.contains(action))
        .map(|action| action.name())
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
    Ok(Parsed::Run(Box::new(args)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(line: &[&str]) -> Vec<String> {
        line.iter().map(|word| word.to_string()).collect()
    }

    /// The options of a line that reads, with the three required ones before `more`.
    fn read(more: &[&str]) -> Args {
        let mut line = vec!["--reference", "r", "--engine", "e=p", "--report", "o"];
        line.extend_from_slice(more);
        match parse(&words(&line)) {
            Ok(Parsed::Run(args)) => *args,
            other => panic!("{line:?} read as {other:?}"),
        }
    }

    fn refused(line: &[&str]) -> String {
        match parse(&words(line)) {
            Err(message) => message,
            other => panic!("{line:?} read as {other:?}"),
        }
    }

    #[test]
    fn every_option_has_its_default() {
        assert_eq!(
            read(&[]),
            Args {
                reference: "r".to_string(),
                engine: vec!["e=p".to_string()],
                report: "o".to_string(),
                sort_memory: "2G".to_string(),
                ..Args::default()
            }
        );
    }

    #[test]
    fn every_option_is_taken_with_its_type() {
        let args = read(&[
            "--engine",
            "f=q",
            "--engine-binary",
            "e=1",
            "--engine-commit",
            "e=2",
            "--plugins",
            "A,B",
            "--run-command",
            "vep --offline",
            "--custom-source",
            "C=c",
            "--custom-source",
            "D=d",
            "--work-dir",
            "w//x/",
            "--residuals",
            "./res",
            "--class-table",
            "t.tsv",
            "--gwas-catalog",
            "g",
            "--revel-table",
            "v",
            "--revel-no-match",
            "--utrannotator-sweep",
            "s1",
            "--utrannotator-sweep",
            "s2/",
            "--json-cache",
            "c/",
            "--spliceai-table",
            "a.synthetic.vcf.gz",
            "--spliceai-table",
            "b.synthetic.vcf.gz",
            "--input",
            "in.vcf",
            "--classes",
            "",
            "--class-population",
            "k=1",
            "--pairing-tool",
            "bin",
            "--pairing-threads",
            " 0_8 ",
            "--sort-memory",
            "16M",
            "--sort-threads",
            "٣",
        ]);
        assert_eq!(args.engine, ["e=p", "f=q"]);
        assert_eq!(
            (args.engine_binary, args.engine_commit),
            (vec!["e=1".to_string()], vec!["e=2".to_string()])
        );
        assert_eq!(args.plugins, "A,B");
        assert_eq!(args.run_command.as_deref(), Some("vep --offline"));
        assert_eq!(args.custom_source, ["C=c", "D=d"]);
        assert_eq!(args.work_dir.as_deref(), Some("w/x"));
        assert_eq!(args.residuals.as_deref(), Some("res"));
        assert_eq!(args.class_table.as_deref(), Some("t.tsv"));
        assert_eq!(
            (args.gwas_catalog.as_deref(), args.revel_table.as_deref()),
            (Some("g"), Some("v"))
        );
        assert!(args.revel_no_match);
        assert_eq!(args.utrannotator_sweep, ["s1", "s2"]);
        assert_eq!(args.json_cache.as_deref(), Some("c"));
        assert_eq!(
            args.spliceai_table,
            ["a.synthetic.vcf.gz", "b.synthetic.vcf.gz"]
        );
        assert_eq!(args.input.as_deref(), Some("in.vcf"));
        assert_eq!(args.classes.as_deref(), Some(""));
        assert_eq!(args.class_population, ["k=1"]);
        assert_eq!(args.pairing_tool.as_deref(), Some("bin"));
        assert_eq!(args.pairing_threads, Some(Integer("8".to_string())));
        assert_eq!(args.sort_memory, "16M");
        assert_eq!(args.sort_threads, Some(Integer("3".to_string())));
    }

    #[test]
    fn option_takes_its_value_after_an_equals_sign_and_by_a_unique_prefix() {
        let args = read(&[
            "--plugins=A",
            "--plug",
            "B",
            "--run=vep --x=1",
            "--class-t=t",
            "--revel-n",
        ]);
        assert_eq!(args.plugins, "B");
        assert_eq!(args.run_command.as_deref(), Some("vep --x=1"));
        assert_eq!(args.class_table.as_deref(), Some("t"));
        assert!(args.revel_no_match);
        assert_eq!(read(&["--classes="]).classes.as_deref(), Some(""));
    }

    #[test]
    fn option_given_again_keeps_its_last_value() {
        let args = read(&[
            "--plugins",
            "A",
            "--plugins",
            "B",
            "--reference",
            "s",
            "--report=p",
        ]);
        assert_eq!(
            (
                args.plugins.as_str(),
                args.reference.as_str(),
                args.report.as_str()
            ),
            ("B", "s", "p")
        );
    }

    #[test]
    fn value_that_opens_with_a_dash_is_taken_with_a_space_or_as_a_negative_number() {
        assert_eq!(
            read(&["--run-command", "--offline --fork 4"])
                .run_command
                .as_deref(),
            Some("--offline --fork 4")
        );
        assert_eq!(
            read(&["--pairing-threads", "-5"]).pairing_threads,
            Some(Integer("-5".to_string()))
        );
        assert_eq!(read(&["--plugins", "-.5x"]).plugins, "-.5x");
        assert_eq!(read(&["--plugins", "-"]).plugins, "-");
        assert_eq!(
            refused(&[
                "--reference",
                "r",
                "--engine",
                "e=p",
                "--report",
                "o",
                "--run-command",
                "--offline"
            ]),
            "argument --run-command: expected one argument"
        );
    }

    #[test]
    fn missing_required_options_are_listed_in_their_order() {
        assert_eq!(
            refused(&[]),
            "the following arguments are required: --reference, --engine, --report"
        );
        assert_eq!(
            refused(&["--report", "o", "--engine", "e=p"]),
            "the following arguments are required: --reference"
        );
        assert_eq!(
            refused(&["--reference", "r", "--bogus", "stray"]),
            "the following arguments are required: --engine, --report"
        );
    }

    #[test]
    fn option_without_its_value_is_refused() {
        assert_eq!(
            refused(&["--reference"]),
            "argument --reference: expected one argument"
        );
        assert_eq!(
            refused(&["--reference", "--engine", "e=p", "--report", "o"]),
            "argument --reference: expected one argument"
        );
        assert_eq!(
            refused(&["--engine", "--", "x"]),
            "argument --engine: expected one argument"
        );
    }

    #[test]
    fn value_that_is_no_whole_number_is_refused() {
        let line = |value: &'static str| {
            vec![
                "--reference",
                "r",
                "--engine",
                "e=p",
                "--report",
                "o",
                "--pairing-threads",
                value,
            ]
        };
        assert_eq!(
            refused(&line("x")),
            "argument --pairing-threads: invalid int value: 'x'"
        );
        assert_eq!(
            refused(&line("1.5")),
            "argument --pairing-threads: invalid int value: '1.5'"
        );
        assert_eq!(
            refused(&line("")),
            "argument --pairing-threads: invalid int value: ''"
        );
        assert_eq!(
            refused(&line("it's")),
            "argument --pairing-threads: invalid int value: \"it's\""
        );
        assert_eq!(
            refused(&["--sort-threads=1__0"]),
            "argument --sort-threads: invalid int value: '1__0'"
        );
        // `int` drops no U+001F around a number, as `str.strip()` does.
        assert_eq!(
            refused(&line("1\u{1f}")),
            "argument --pairing-threads: invalid int value: '1\\x1f'"
        );
    }

    #[test]
    fn arguments_no_option_takes_are_listed() {
        let base = ["--reference", "r", "--engine", "e=p", "--report", "o"];
        let with = |more: &[&'static str]| refused(&[&base[..], more].concat());
        assert_eq!(with(&["--bogus"]), "unrecognized arguments: --bogus");
        assert_eq!(
            with(&["--bogus=1", "stray", "--other"]),
            "unrecognized arguments: --bogus=1 stray --other"
        );
        assert_eq!(with(&["-x"]), "unrecognized arguments: -x");
        assert_eq!(with(&["-"]), "unrecognized arguments: -");
        assert_eq!(
            with(&["--", "--plugins", "x"]),
            "unrecognized arguments: -- --plugins x"
        );
        assert_eq!(
            refused(&[&["stray"][..], &base[..], &["more"][..]].concat()),
            "unrecognized arguments: stray more"
        );
    }

    #[test]
    fn prefix_of_several_options_names_them_all() {
        assert_eq!(
            refused(&["--re", "x"]),
            "ambiguous option: --re could match --reference, --report, --residuals, --revel-table, --revel-no-match"
        );
        assert_eq!(
            refused(&["--engine-=x"]),
            "ambiguous option: --engine-=x could match --engine-binary, --engine-commit"
        );
        assert_eq!(
            refused(&["--sort"]),
            "ambiguous option: --sort could match --sort-memory, --sort-threads"
        );
    }

    #[test]
    fn word_that_opens_with_a_dash_and_an_equals_sign_fits_every_option() {
        const EVERY_OPTION: &str = "-h, --help, --reference, --engine, --engine-binary, --engine-commit, --plugins, --run-command, --custom-source, --work-dir, --report, --residuals, --class-table, --gwas-catalog, --revel-table, --revel-no-match, --utrannotator-sweep, --json-cache, --spliceai-table, --input, --classes, --class-population, --pairing-tool, --pairing-threads, --sort-memory, --sort-threads";
        for word in ["-=x", "-=", "-=x y", "-=-1"] {
            assert_eq!(
                refused(&["--plugins", "A", word]),
                format!("ambiguous option: {word} could match {EVERY_OPTION}")
            );
        }
        assert_eq!(
            refused(&["--reference", "-=x"]),
            "argument --reference: expected one argument"
        );
        assert_eq!(parse(&words(&["-h", "-=x"])), Ok(Parsed::Help));
    }

    #[test]
    fn option_that_takes_no_value_refuses_one() {
        assert_eq!(
            refused(&["--revel-no-match=1"]),
            "argument --revel-no-match: ignored explicit argument '1'"
        );
        assert_eq!(
            refused(&["--help=1"]),
            "argument -h/--help: ignored explicit argument '1'"
        );
        assert_eq!(
            refused(&["-h=1"]),
            "argument -h/--help: ignored explicit argument '1'"
        );
        assert_eq!(
            refused(&["--help="]),
            "argument -h/--help: ignored explicit argument ''"
        );
    }

    #[test]
    fn help_is_taken_where_it_stands() {
        for line in [
            &["-h"][..],
            &["--help"],
            &["--he"],
            &["-hh"],
            &["-hx"],
            &["--help", "--pairing-threads", "x"],
            &["--bogus", "-h"],
        ] {
            assert_eq!(parse(&words(line)), Ok(Parsed::Help), "{line:?}");
        }
        assert_eq!(
            refused(&["--pairing-threads", "x", "--help"]),
            "argument --pairing-threads: invalid int value: 'x'"
        );
        // After an option that takes a value, the help's name is no value.
        assert_eq!(
            refused(&["--reference", "-h"]),
            "argument --reference: expected one argument"
        );
    }
}
