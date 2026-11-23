// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `compare` command line, read as CPython's `argparse` reads a parser of these
//! options: `--name value` and `--name=value`, a unique prefix of a long option, a value
//! that opens with `-` only where it reads as a negative number, and `--` before
//! arguments that are no options.

use super::text::{decimal_digit, push_str_repr, py_float, py_int, PyInt};

pub const USAGE: &str = "usage: vep-plugin-pairs compare [-h] [--perl-output PATH] [--rust-output PATH] [--perl-dir DIR] [--rust-dir DIR] --plugins NAMES [--fields NAMES] [--json-output PATH] [--max-mismatches N] [--tolerance X] [--min-annotated N] [--threads N] [--block-bytes N] [--shard-bytes N] [--work-dir DIR]";
pub const HELP: &str = "\
Compares the fields plugins write between an Ensembl VEP output and a vep-rs output.

options:
  -h, --help          show this help and exit
  --perl-output PATH  one Ensembl VEP output, default format or VCF, plain or gzip
  --rust-output PATH  one vep-rs output of the same input
  --perl-dir DIR      a directory of Ensembl VEP outputs
  --rust-dir DIR      a directory of vep-rs outputs; a file pairs with the one of the same
                      name without its ending and a leading perl_ or rust_
  --plugins NAMES     comma-separated plugin names, or all for every plugin of the registry
  --fields NAMES      comma-separated field names compared in place of the registry's
                      entries, for a run that names one plugin
  --json-output PATH  write the results as JSON to this file
  --max-mismatches N  mismatches kept in the JSON output per plugin and file pair
                      (default 10; 0 keeps none); the printed report shows three at most
  --tolerance X       two numbers less than this apart are equal (default 1e-06)
  --min-annotated N   records each engine must annotate, or the plugin fails (default 1):
                      a plugin that annotates nothing has not been measured
  --threads N         threads of the run (default: every core)
  --block-bytes N     bytes of a plain file a thread reads at a time (default 67108864)
  --shard-bytes N     bytes of the longer output a shard of the scratch is sized for
                      (default 134217728)
  --work-dir DIR      directory the run's scratch is made in (default: the system's
                      temporary directory)

VEP_PLUGIN_PAIRS_TIMINGS in the environment prints the time of each step of a file pair
on standard error.";

/// The options `argparse` knows, in the order it lists them where a prefix fits several.
const KNOWN: [&str; 12] = [
    "-h",
    "--help",
    "--perl-output",
    "--rust-output",
    "--perl-dir",
    "--rust-dir",
    "--plugins",
    "--fields",
    "--json-output",
    "--max-mismatches",
    "--tolerance",
    "--min-annotated",
];
/// The tool's own options, taken by their whole names alone so that every prefix means what
/// it means to `argparse`.
const THREADS: &str = "--threads";
const BLOCK_BYTES: &str = "--block-bytes";
const SHARD_BYTES: &str = "--shard-bytes";
const WORK_DIR: &str = "--work-dir";

#[derive(Debug, PartialEq)]
pub struct Options {
    pub perl_output: Option<String>,
    pub rust_output: Option<String>,
    pub perl_dir: Option<String>,
    pub rust_dir: Option<String>,
    pub plugins: String,
    pub fields: Option<String>,
    pub json_output: Option<String>,
    /// At most this many mismatches are kept; a count beyond `usize` keeps every one.
    pub max_mismatches: usize,
    pub tolerance: f64,
    pub min_annotated: PyInt,
    pub threads: usize,
    pub block_bytes: u64,
    pub shard_bytes: u64,
    pub work_dir: Option<String>,
}

#[derive(Debug, PartialEq)]
pub enum Command {
    Run(Box<Options>),
    Help,
    /// A wrong command line, with `argparse`'s message for it.
    Refused(String),
}

enum Word {
    /// No option: a value, or an argument no option takes.
    Argument,
    /// `--`, after which every word is an argument.
    Separator,
    /// An option, with each known option its text fits and the value written with it.
    Option(Vec<(&'static str, Option<String>)>),
    /// An option the parser does not know.
    Unknown,
}

fn known(name: &str) -> Option<&'static str> {
    KNOWN
        .iter()
        .copied()
        .chain([THREADS, BLOCK_BYTES, SHARD_BYTES, WORK_DIR])
        .find(|option| *option == name)
}

fn classify(word: &str) -> Word {
    if !word.starts_with('-') {
        return Word::Argument;
    }
    if let Some(option) = known(word) {
        return Word::Option(vec![(option, None)]);
    }
    if word.len() == 1 {
        return Word::Argument;
    }
    let (name, explicit) = match word.split_once('=') {
        Some((name, explicit)) => (name, Some(explicit.to_string())),
        None => (word, None),
    };
    if explicit.is_some() {
        if let Some(option) = known(name) {
            return Word::Option(vec![(option, explicit)]);
        }
    }
    // Behind a single dash the text before `=` is a prefix too, and the dash alone, the
    // word `-=value`, is the prefix of every option.
    let fitting: Vec<(&'static str, Option<String>)> = if word.starts_with("--") || name == "-" {
        KNOWN
            .iter()
            .filter(|option| option.starts_with(name))
            .map(|option| (*option, explicit.clone()))
            .collect()
    } else if let Some(rest) = word.strip_prefix("-h") {
        // A single dash and a letter is the short option with the rest as its value.
        vec![("-h", Some(rest.to_string()))]
    } else {
        Vec::new()
    };
    if !fitting.is_empty() {
        return Word::Option(fitting);
    }
    let mut after_dash = word[1..].chars();
    let negative_number = match after_dash.next() {
        Some('.') => after_dash
            .next()
            .is_some_and(|c| decimal_digit(c).is_some()),
        Some(c) => decimal_digit(c).is_some(),
        None => false,
    };
    if negative_number || word.contains(' ') {
        return Word::Argument;
    }
    Word::Unknown
}

fn repr(text: &str) -> String {
    let mut out = String::new();
    push_str_repr(&mut out, text);
    out
}

/// What `argparse` makes of text written against the help option. Behind `=`, and on the
/// long option, the text is a value the option does not take. Letters behind the short
/// option are further short options, so the help is given, unless a dash follows or an `=`
/// follows a repeated `h`.
fn help_with_text(word: &str, text: &str) -> Result<Command, String> {
    let ignored = |text: &str| {
        Err(format!(
            "argument -h/--help: ignored explicit argument {}",
            repr(text)
        ))
    };
    if word.starts_with("--") || word.starts_with("-h=") {
        return ignored(text);
    }
    let mut rest = text;
    loop {
        if rest.starts_with('-') {
            return ignored(rest);
        }
        match rest.strip_prefix('h') {
            Some("") | None => return Ok(Command::Help),
            Some(after) => match after.strip_prefix('=') {
                Some(value) => return ignored(value),
                None => rest = after,
            },
        }
    }
}

/// Reads the words after `compare`.
pub fn parse(words: &[String]) -> Command {
    match read(words) {
        Ok(command) => command,
        Err(message) => Command::Refused(message),
    }
}

fn read(words: &[String]) -> Result<Command, String> {
    let mut kinds: Vec<Word> = Vec::with_capacity(words.len());
    let mut separated = false;
    for word in words {
        kinds.push(if separated {
            Word::Argument
        } else if word == "--" {
            separated = true;
            Word::Separator
        } else {
            classify(word)
        });
    }

    let mut perl_output = None;
    let mut rust_output = None;
    let mut perl_dir = None;
    let mut rust_dir = None;
    let mut plugins = None;
    let mut fields = None;
    let mut json_output = None;
    let mut max_mismatches = PyInt::from_i64(10);
    let mut tolerance = 1e-6;
    let mut min_annotated = PyInt::from_i64(1);
    let mut threads = std::thread::available_parallelism().map_or(1, usize::from);
    let mut block_bytes = 1 << 26;
    let mut shard_bytes = 1 << 27;
    let mut work_dir = None;
    let mut unrecognized: Vec<&str> = Vec::new();

    let mut at = 0;
    while at < words.len() {
        let fitting = match &kinds[at] {
            Word::Argument | Word::Separator | Word::Unknown => {
                unrecognized.push(&words[at]);
                at += 1;
                continue;
            }
            Word::Option(fitting) => fitting,
        };
        if fitting.len() > 1 {
            let names: Vec<&str> = fitting.iter().map(|(option, _)| *option).collect();
            return Err(format!(
                "ambiguous option: {} could match {}",
                words[at],
                names.join(", ")
            ));
        }
        let (option, explicit) = &fitting[0];
        if matches!(*option, "-h" | "--help") {
            return match explicit.as_deref() {
                None => Ok(Command::Help),
                Some(text) => help_with_text(&words[at], text),
            };
        }
        let value: &str = match explicit {
            Some(value) => {
                at += 1;
                value
            }
            None => {
                if !matches!(kinds.get(at + 1), Some(Word::Argument)) {
                    return Err(format!("argument {option}: expected one argument"));
                }
                at += 2;
                &words[at - 1]
            }
        };
        let invalid =
            |kind: &str| format!("argument {option}: invalid {kind} value: {}", repr(value));
        match *option {
            "--perl-output" => perl_output = Some(value.to_string()),
            "--rust-output" => rust_output = Some(value.to_string()),
            "--perl-dir" => perl_dir = Some(value.to_string()),
            "--rust-dir" => rust_dir = Some(value.to_string()),
            "--plugins" => plugins = Some(value.to_string()),
            "--fields" => fields = Some(value.to_string()),
            "--json-output" => json_output = Some(value.to_string()),
            "--max-mismatches" => max_mismatches = py_int(value).ok_or_else(|| invalid("int"))?,
            "--tolerance" => tolerance = py_float(value).ok_or_else(|| invalid("float"))?,
            "--min-annotated" => min_annotated = py_int(value).ok_or_else(|| invalid("int"))?,
            WORK_DIR => work_dir = Some(value.to_string()),
            _ => {
                let count = py_int(value).ok_or_else(|| invalid("int"))?;
                let positive = u64::try_from(count.saturating_i128())
                    .ok()
                    .filter(|count| *count >= 1)
                    .ok_or_else(|| {
                        format!("argument {option}: must be 1 or greater, got {count}")
                    })?;
                match *option {
                    THREADS => threads = usize::try_from(positive).unwrap_or(usize::MAX),
                    BLOCK_BYTES => block_bytes = positive,
                    _ => shard_bytes = positive,
                }
            }
        }
    }
    let Some(plugins) = plugins else {
        return Err("the following arguments are required: --plugins".to_string());
    };
    if !unrecognized.is_empty() {
        return Err(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }
    if max_mismatches.is_negative() {
        // As the bound of the kept mismatches a negative value would drop them from the end.
        return Err(format!(
            "argument --max-mismatches: must be 0 or greater, got {max_mismatches}"
        ));
    }
    Ok(Command::Run(Box::new(Options {
        perl_output,
        rust_output,
        perl_dir,
        rust_dir,
        plugins,
        fields,
        json_output,
        max_mismatches: usize::try_from(max_mismatches.saturating_i128()).unwrap_or(usize::MAX),
        tolerance,
        min_annotated,
        threads,
        block_bytes,
        shard_bytes,
        work_dir,
    })))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(words: &[&str]) -> Command {
        let words: Vec<String> = words.iter().map(|word| word.to_string()).collect();
        parse(&words)
    }

    fn options(words: &[&str]) -> Options {
        match parsed(words) {
            Command::Run(options) => *options,
            other => panic!("{words:?} gave {other:?}"),
        }
    }

    fn refusal(words: &[&str]) -> String {
        match parsed(words) {
            Command::Refused(message) => message,
            other => panic!("{words:?} gave {other:?}"),
        }
    }

    #[test]
    fn every_option_has_its_default() {
        let read = options(&["--plugins", "REVEL"]);
        assert_eq!(read.plugins, "REVEL");
        assert_eq!(
            (
                &read.perl_output,
                &read.rust_output,
                &read.perl_dir,
                &read.rust_dir
            ),
            (&None, &None, &None, &None)
        );
        assert_eq!((&read.fields, &read.json_output), (&None, &None));
        assert_eq!(read.max_mismatches, 10);
        assert_eq!(read.tolerance, 1e-6);
        assert_eq!(read.min_annotated, PyInt::from_i64(1));
        assert!(read.threads >= 1);
        assert_eq!(read.block_bytes, 1 << 26);
        assert_eq!(read.shard_bytes, 1 << 27);
        assert_eq!(read.work_dir, None);
    }

    #[test]
    fn options_are_read_with_a_separate_value_an_equals_sign_or_a_unique_prefix() {
        let read = options(&[
            "--perl-output",
            "p.txt",
            "--rust-o=r.txt",
            "--plug",
            "CADD,REVEL",
            "--f=a,b",
            "--json",
            "out.json",
            "--max",
            "3",
            "--t",
            "1e-4",
            "--min=0",
            "--threads",
            "2",
            "--block-bytes=4096",
            "--shard-bytes",
            "64",
            "--work-dir=scratch here",
        ]);
        assert_eq!(read.perl_output.as_deref(), Some("p.txt"));
        assert_eq!(read.rust_output.as_deref(), Some("r.txt"));
        assert_eq!(read.plugins, "CADD,REVEL");
        assert_eq!(read.fields.as_deref(), Some("a,b"));
        assert_eq!(read.json_output.as_deref(), Some("out.json"));
        assert_eq!(read.max_mismatches, 3);
        assert_eq!(read.tolerance, 1e-4);
        assert_eq!(read.min_annotated, PyInt::from_i64(0));
        assert_eq!(read.threads, 2);
        assert_eq!(read.block_bytes, 4096);
        assert_eq!(read.shard_bytes, 64);
        assert_eq!(read.work_dir.as_deref(), Some("scratch here"));
        let read = options(&[
            "--plugins",
            "a",
            "--plugins=b",
            "--min-annotated",
            "-5",
            "--tolerance",
            "-.5",
        ]);
        assert_eq!(read.plugins, "b");
        assert_eq!(read.min_annotated, PyInt::from_i64(-5));
        assert_eq!(read.tolerance, -0.5);
        assert_eq!(options(&["--plugins="]).plugins, "");
        assert_eq!(options(&["--plugins", "a b"]).plugins, "a b");
        assert_eq!(options(&["--plugins", "-"]).plugins, "-");
        assert_eq!(
            options(&["--plugins", "x", "--max-mismatches", " 1_0 "]).max_mismatches,
            10
        );
    }

    #[test]
    fn a_wrong_command_line_is_refused_with_the_message_argparse_gives() {
        let cases: [(&[&str], &str); 21] = [
            (&[], "the following arguments are required: --plugins"),
            (&["--plugins"], "argument --plugins: expected one argument"),
            (
                &["--plugins", "--fields", "x"],
                "argument --plugins: expected one argument",
            ),
            (
                &["--plugins", "-x"],
                "argument --plugins: expected one argument",
            ),
            (
                &["--plugins", "--"],
                "argument --plugins: expected one argument",
            ),
            (
                &["--perl", "x", "--plugins", "y"],
                "ambiguous option: --perl could match --perl-output, --perl-dir",
            ),
            (
                &["--plugins", "x", "--m=1"],
                "ambiguous option: --m=1 could match --max-mismatches, --min-annotated",
            ),
            (
                &["--plugins", "x", "extra", "--bogus", "more"],
                "unrecognized arguments: extra --bogus more",
            ),
            (
                &["--plugins", "x", "--thr", "2"],
                "unrecognized arguments: --thr 2",
            ),
            (
                &["--plugins", "x", "--max-mismatches", "many"],
                "argument --max-mismatches: invalid int value: 'many'",
            ),
            (
                &["--plugins", "x", "--max-mismatches", "1.5"],
                "argument --max-mismatches: invalid int value: '1.5'",
            ),
            (
                &["--plugins", "x", "--tolerance", "it's"],
                "argument --tolerance: invalid float value: \"it's\"",
            ),
            (
                &["--plugins", "x", "--max-mismatches", "-1"],
                "argument --max-mismatches: must be 0 or greater, got -1",
            ),
            (
                &["--plugins", "x", "--max-mismatches=-007"],
                "argument --max-mismatches: must be 0 or greater, got -7",
            ),
            (
                &["--plugins", "x", "--help=1"],
                "argument -h/--help: ignored explicit argument '1'",
            ),
            (
                &["--plugins", "x", "--threads", "0"],
                "argument --threads: must be 1 or greater, got 0",
            ),
            (
                &["--plugins", "x", "--block-bytes", "-4"],
                "argument --block-bytes: must be 1 or greater, got -4",
            ),
            (
                &["--plugins", "x", "--block", "4"],
                "unrecognized arguments: --block 4",
            ),
            (
                &["--plugins", "x", "--shard-bytes", "0"],
                "argument --shard-bytes: must be 1 or greater, got 0",
            ),
            (
                &["--plugins", "x", "--work-dir"],
                "argument --work-dir: expected one argument",
            ),
            (
                &["--plugins", "x", "--work", "here"],
                "unrecognized arguments: --work here",
            ),
        ];
        for (words, message) in cases {
            assert_eq!(refusal(words), message, "{words:?}");
        }
    }

    #[test]
    fn a_word_that_opens_with_a_dash_and_an_equals_sign_fits_every_option() {
        const EVERY_OPTION: &str = "-h, --help, --perl-output, --rust-output, --perl-dir, --rust-dir, --plugins, --fields, --json-output, --max-mismatches, --tolerance, --min-annotated";
        for (words, word) in [
            (&["--plugins", "x", "-=x"][..], "-=x"),
            (&["-=x"], "-=x"),
            (&["-="], "-="),
            (&["--plugins", "x", "-=x y"], "-=x y"),
            (&["--plugins", "x", "-=-1"], "-=-1"),
            (&["extra", "-=h", "--plugins", "x"], "-=h"),
            (&["--plugins", "x", "-=--plugins", "-h"], "-=--plugins"),
        ] {
            assert_eq!(
                refusal(words),
                format!("ambiguous option: {word} could match {EVERY_OPTION}"),
                "{words:?}"
            );
        }
        // Such a word is an option, so no option takes it as its value.
        assert_eq!(
            refusal(&["--plugins", "-=x"]),
            "argument --plugins: expected one argument"
        );
        assert_eq!(parsed(&["-h", "-=x"]), Command::Help);
        // Behind `--` it is an argument, and with other text before the `=` it fits nothing.
        assert_eq!(
            refusal(&["--plugins", "x", "--", "-=x"]),
            "unrecognized arguments: -- -=x"
        );
        assert_eq!(
            refusal(&["--plugins", "x", "-x=1", "--=x"]),
            "ambiguous option: --=x could match --help, --perl-output, --rust-output, --perl-dir, --rust-dir, --plugins, --fields, --json-output, --max-mismatches, --tolerance, --min-annotated"
        );
        assert_eq!(
            refusal(&["--plugins", "x", "-x=1"]),
            "unrecognized arguments: -x=1"
        );
    }

    #[test]
    fn help_is_asked_by_either_name_and_stops_the_reading() {
        assert_eq!(parsed(&["-h"]), Command::Help);
        assert_eq!(parsed(&["--help", "--bogus"]), Command::Help);
        assert_eq!(parsed(&["--he"]), Command::Help);
        assert_eq!(
            parsed(&["--plugins", "x", "-h", "--max-mismatches", "many"]),
            Command::Help
        );
        assert!(matches!(
            parsed(&["--max-mismatches", "many", "-h"]),
            Command::Refused(_)
        ));
    }

    #[test]
    fn negative_max_mismatches_is_refused_by_name() {
        let message = refusal(&[
            "--perl-output",
            "p",
            "--rust-output",
            "r",
            "--plugins",
            "REVEL",
            "--json-output",
            "o.json",
            "--max-mismatches",
            "-1",
        ]);
        assert!(message.contains("argument --max-mismatches"));
    }
}
