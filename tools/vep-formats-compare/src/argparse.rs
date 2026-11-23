// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line's options read as `argparse` of the interpreter reads them.
//!
//! The options are declared as that `argparse` takes them (`cli::build_parser`). This
//! module reads a command line as that parser does: which argument is an option, a value
//! or neither; an option spelt by an unambiguous prefix; the `--option=value` form; the
//! order in which the errors are met, and their texts. The help and usage texts are that
//! parser's own at 80 columns, stored with the tool: the width and the colours of a
//! terminal are not read.
//!
//! An argument that is not UTF-8 is held as the interpreter holds it, each byte that is no
//! part of a character as one code point of its own (`text_of`), so that the file it names
//! is opened by its bytes and every text that quotes it prints the byte as the interpreter
//! does (`standard_output`, `standard_error`, `repr`, and the report's writer).

use crate::pyfmt;
use std::collections::BTreeSet;
use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

/// The code points that stand for the bytes 0x80 to 0xFF of an argument that is not UTF-8.
/// The interpreter holds such a byte as the lone surrogate U+DC80 to U+DCFF, which a
/// `String` cannot hold; these are the code points of the last plane that end in the same
/// sixteen bits.
const HELD_BYTES: std::ops::RangeInclusive<u32> = 0x10DC80..=0x10DCFF;

/// Whether `text_of` has held a byte in this run. Without one no text is read as holding
/// bytes, and a character of `HELD_BYTES` in an output is printed as the character it is.
static HOLDS_BYTES: AtomicBool = AtomicBool::new(false);

fn held(byte: u8) -> char {
    char::from_u32(0x10DC00 + u32::from(byte)).unwrap_or(char::REPLACEMENT_CHARACTER)
}

fn held_byte(c: char) -> Option<u8> {
    HELD_BYTES
        .contains(&(c as u32))
        .then(|| (c as u32 - 0x10DC00) as u8)
}

/// The text of an argument: its characters, and one code point of `HELD_BYTES` for each
/// byte that is no part of a UTF-8 character. A character of that range that the argument
/// spells in UTF-8 is held as its four bytes, so that `bytes_of` returns every argument's
/// own bytes.
pub fn text_of(arg: &OsStr) -> String {
    let mut bytes = arg.as_bytes();
    let mut text = String::with_capacity(bytes.len());
    let mut holds = false;
    while !bytes.is_empty() {
        let (valid, bad) = match std::str::from_utf8(bytes) {
            Ok(valid) => (valid, 0),
            Err(error) => {
                let valid = std::str::from_utf8(&bytes[..error.valid_up_to()]).unwrap_or("");
                let bad = error
                    .error_len()
                    .unwrap_or(bytes.len() - error.valid_up_to());
                (valid, bad)
            }
        };
        for c in valid.chars() {
            if held_byte(c).is_some() {
                holds = true;
                text.extend(c.encode_utf8(&mut [0; 4]).bytes().map(held));
            } else {
                text.push(c);
            }
        }
        let rest = &bytes[valid.len()..];
        holds |= bad > 0;
        text.extend(rest[..bad].iter().copied().map(held));
        bytes = &rest[bad..];
    }
    if holds {
        HOLDS_BYTES.store(true, Ordering::Relaxed);
    }
    text
}

/// Whether an argument of this run held a byte that is no part of a UTF-8 character.
pub fn holds_bytes() -> bool {
    HOLDS_BYTES.load(Ordering::Relaxed)
}

/// The bytes a text of `text_of` stands for.
pub fn bytes_of(text: &str) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(text.len());
    for c in text.chars() {
        match held_byte(c) {
            Some(byte) => bytes.push(byte),
            None => bytes.extend_from_slice(c.encode_utf8(&mut [0; 4]).as_bytes()),
        }
    }
    bytes
}

/// A text of `text_of` with each held byte `0xNN` spelt as the six characters `\udcNN`.
pub fn backslashreplace(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match held_byte(c) {
            Some(byte) => out.push_str(&format!("\\udc{byte:02x}")),
            None => out.push(c),
        }
    }
    out
}

/// What the interpreter writes of `text` to the standard output, whose encoder turns a
/// held byte back into itself.
pub fn standard_output(text: &str) -> Vec<u8> {
    if holds_bytes() {
        bytes_of(text)
    } else {
        text.as_bytes().to_vec()
    }
}

/// What the interpreter writes of `text` to the standard error, whose encoder spells a
/// held byte `0xNN` as `\udcNN`.
pub fn standard_error(text: &str) -> String {
    if holds_bytes() {
        backslashreplace(text)
    } else {
        text.to_string()
    }
}

/// `text` with every escape that opens with `from` and ends in two hexadecimal digits of
/// 0x80 or more rewritten to open with `to`. `text` is a `repr()` or a JSON text, where a
/// backslash opens an escape and the character after it belongs to that escape.
pub fn respell_escapes(text: &str, from: &str, to: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find('\\') {
        out.push_str(&rest[..at]);
        let escape = &rest[at..];
        let digits = escape
            .strip_prefix(from)
            .and_then(|digits| digits.get(..2))
            .filter(|digits| u8::from_str_radix(digits, 16).is_ok_and(|byte| byte >= 0x80));
        match digits {
            Some(digits) => {
                out.push_str(to);
                out.push_str(digits);
                rest = &escape[from.len() + 2..];
            }
            None => {
                let mut chars = escape.chars();
                out.extend(chars.by_ref().take(2));
                rest = chars.as_str();
            }
        }
    }
    out.push_str(rest);
    out
}

/// `repr(text)` of a text of `text_of`: a held byte `0xNN` is spelt `\udcNN`.
pub fn repr(text: &str) -> String {
    let plain = pyfmt::repr_str(text);
    if text.chars().any(|c| held_byte(c).is_some()) {
        respell_escapes(&plain, "\\U0010dc", "\\udc")
    } else {
        plain
    }
}

/// The root and the parts of `pathlib.Path(text)`: one `/` or exactly two as the root,
/// and the parts between `/` that are not empty and not `.`.
fn path_parts(raw: &[u8]) -> (&'static [u8], Vec<&[u8]>) {
    let root: &[u8] = match raw.iter().take_while(|byte| **byte == b'/').count() {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let parts = raw
        .split(|byte| *byte == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
        .collect();
    (root, parts)
}

fn path_from_parts(root: &[u8], parts: &[&[u8]]) -> PathBuf {
    let mut raw = root.to_vec();
    raw.extend_from_slice(&parts.join(&b'/'));
    if raw.is_empty() {
        raw.push(b'.');
    }
    PathBuf::from(OsString::from_vec(raw))
}

/// `pathlib.Path(text)` of an argument's text, as the path `str()` prints and every
/// reader opens: repeated `/`, `.` parts and a trailing `/` dropped, `..` kept, exactly two
/// leading `/` kept, the empty text as `.`.
pub fn path_of(text: &str) -> PathBuf {
    let raw = bytes_of(text);
    let (root, parts) = path_parts(&raw);
    path_from_parts(root, &parts)
}

/// `str(path)` of a path of `path_of`, as a text of `text_of`.
pub fn path_text(path: &Path) -> String {
    text_of(path.as_os_str())
}

/// `path.name` of a path of `path_of`: its last part, the empty text where it has none
/// (`.`, `/`); `..` is a part.
pub fn path_name(path: &Path) -> String {
    let (_, parts) = path_parts(path.as_os_str().as_bytes());
    parts
        .last()
        .map(|name| text_of(OsStr::from_bytes(name)))
        .unwrap_or_default()
}

/// `path.parent` of a path of `path_of`: the path without its last part, itself where it
/// has none.
pub fn path_parent(path: &Path) -> PathBuf {
    let (root, mut parts) = path_parts(path.as_os_str().as_bytes());
    parts.pop();
    path_from_parts(root, &parts)
}

/// What an option does with the command line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// Prints the help and ends the run.
    Help,
    /// A switch: it takes no value.
    StoreTrue,
    /// Takes one value; a later occurrence replaces an earlier one.
    Store,
    /// Takes one value an occurrence; the values are kept in order.
    Append,
}

impl Action {
    fn takes_a_value(self) -> bool {
        matches!(self, Action::Store | Action::Append)
    }
}

/// What the parser makes of an option's value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// The text as given.
    Text,
    /// `pathlib.Path(text)`.
    Path,
    /// `int(text)`.
    Int,
}

/// One option of the parser.
#[derive(Debug, Clone)]
pub struct Argument {
    pub option_strings: Vec<&'static str>,
    /// The name its value is kept under: the first long option string without its dashes,
    /// `-` as `_`.
    pub dest: String,
    pub action: Action,
    pub kind: Kind,
    pub choices: Option<&'static [&'static str]>,
}

impl Argument {
    /// The name an error gives the option: its option strings joined by `/`.
    fn name(&self) -> String {
        self.option_strings.join("/")
    }
}

/// The value an occurrence of an option was given.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    /// A switch was met.
    Switch,
    Text(String),
    Path(PathBuf),
    /// `int(text)`; a number beyond `i64` is held as the nearest `i64`.
    Int(i64),
}

/// One occurrence of an option on the command line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Taken {
    /// The option's place among the parser's arguments.
    pub argument: usize,
    pub value: Value,
}

/// How the parser ends a run by itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Exit {
    /// The help option was met: the help goes to the standard output, exit code 0.
    Help,
    /// A usage error: the usage and the error line go to the standard error, exit code 2.
    Error(String),
}

/// The command line's options with the texts of their usage and help.
#[derive(Debug)]
pub struct ArgumentParser {
    pub prog: &'static str,
    usage: &'static str,
    help: &'static str,
    arguments: Vec<Argument>,
    /// The groups of options of which a command line may name one.
    exclusive: Vec<Vec<usize>>,
}

/// One reading of an argument as an option: the option, absent for one the parser does not
/// have; its option string; and the value written in the same argument with the separator
/// before it (`=`, or nothing after a one-letter option).
#[derive(Debug, Clone)]
struct OptionTuple {
    argument: Option<usize>,
    option_string: String,
    sep: Option<&'static str>,
    explicit_arg: Option<String>,
}

/// `text.partition("=")`, the separator and the tail absent where the text has no `=`.
fn split_at_equals(text: &str) -> (&str, Option<&'static str>, Option<String>) {
    match text.split_once('=') {
        Some((head, tail)) => (head, Some("="), Some(tail.to_string())),
        None => (text, None, None),
    }
}

/// Whether the text opens as a negative number does (`-\.?\d`): the parser reads such an
/// argument as a value, having no option of that shape.
fn opens_as_a_negative_number(text: &str) -> bool {
    let mut chars = text.chars();
    if chars.next() != Some('-') {
        return false;
    }
    match chars.next() {
        Some('.') => chars.next().is_some_and(pyfmt::char_is_decimal),
        Some(c) => pyfmt::char_is_decimal(c),
        None => false,
    }
}

impl ArgumentParser {
    /// A parser with its help option; `usage` and `help` are the texts `format_usage` and
    /// `format_help` return.
    pub fn new(prog: &'static str, usage: &'static str, help: &'static str) -> Self {
        let mut parser = ArgumentParser {
            prog,
            usage,
            help,
            arguments: Vec::new(),
            exclusive: Vec::new(),
        };
        parser.add(vec!["-h", "--help"], Action::Help, Kind::Text, None);
        parser
    }

    fn add(
        &mut self,
        option_strings: Vec<&'static str>,
        action: Action,
        kind: Kind,
        choices: Option<&'static [&'static str]>,
    ) -> usize {
        let long = option_strings
            .iter()
            .find(|option_string| option_string.starts_with("--"))
            .unwrap_or(&option_strings[0]);
        let dest = long.trim_start_matches('-').replace('-', "_");
        self.arguments.push(Argument {
            option_strings,
            dest,
            action,
            kind,
            choices,
        });
        self.arguments.len() - 1
    }

    /// Adds the option `option_string` and returns its place among the arguments.
    pub fn add_argument(
        &mut self,
        option_string: &'static str,
        action: Action,
        kind: Kind,
        choices: Option<&'static [&'static str]>,
    ) -> usize {
        self.add(vec![option_string], action, kind, choices)
    }

    /// States that a command line may name one of `arguments` alone.
    pub fn add_mutually_exclusive_group(&mut self, arguments: &[usize]) {
        self.exclusive.push(arguments.to_vec());
    }

    pub fn arguments(&self) -> &[Argument] {
        &self.arguments
    }

    /// The name the value of the argument at `argument` is kept under.
    pub fn dest(&self, argument: usize) -> &str {
        &self.arguments[argument].dest
    }

    /// The usage block.
    pub fn format_usage(&self) -> &'static str {
        self.usage
    }

    /// The whole help.
    pub fn format_help(&self) -> &'static str {
        self.help
    }

    /// What a usage error writes to the standard error: the usage block and the error line.
    pub fn error_text(&self, message: &str) -> String {
        format!("{}{}: error: {}\n", self.usage, self.prog, message)
    }

    /// Every option string with its argument, in the order the options were added.
    fn option_strings(&self) -> impl Iterator<Item = (usize, &'static str)> + '_ {
        self.arguments
            .iter()
            .enumerate()
            .flat_map(|(at, argument)| argument.option_strings.iter().map(move |s| (at, *s)))
    }

    fn option(&self, option_string: &str) -> Option<usize> {
        self.option_strings()
            .find(|(_, known)| *known == option_string)
            .map(|(at, _)| at)
    }

    fn argument_error(&self, argument: usize, message: &str) -> Exit {
        Exit::Error(format!(
            "argument {}: {message}",
            self.arguments[argument].name()
        ))
    }

    /// The readings of `arg` as an option, `None` for an argument that is a value: one
    /// that does not open with `-`, `-` alone, and, where no option matches it, one that
    /// opens as a negative number or holds a space.
    fn parse_optional(&self, arg: &str) -> Option<Vec<OptionTuple>> {
        if !arg.starts_with('-') {
            return None;
        }
        if let Some(argument) = self.option(arg) {
            return Some(vec![OptionTuple {
                argument: Some(argument),
                option_string: arg.to_string(),
                sep: None,
                explicit_arg: None,
            }]);
        }
        if arg == "-" {
            return None;
        }
        let (option_string, sep, explicit_arg) = split_at_equals(arg);
        if sep.is_some() {
            if let Some(argument) = self.option(option_string) {
                return Some(vec![OptionTuple {
                    argument: Some(argument),
                    option_string: option_string.to_string(),
                    sep,
                    explicit_arg,
                }]);
            }
        }
        let option_tuples = self.option_tuples(arg);
        if !option_tuples.is_empty() {
            return Some(option_tuples);
        }
        if opens_as_a_negative_number(arg) || arg.contains(' ') {
            return None;
        }
        Some(vec![OptionTuple {
            argument: None,
            option_string: arg.to_string(),
            sep: None,
            explicit_arg: None,
        }])
    }

    /// The options `arg` may abbreviate: after two dashes every option string that opens
    /// with the text before its first `=`; after one dash the one-letter option of its
    /// first two characters, the rest being that option's value, and every option string
    /// that opens with the text before the first `=`.
    fn option_tuples(&self, arg: &str) -> Vec<OptionTuple> {
        let (option_prefix, sep, explicit_arg) = split_at_equals(arg);
        let mut result = Vec::new();
        if arg[1..].starts_with('-') {
            for (argument, option_string) in self.option_strings() {
                if option_string.starts_with(option_prefix) {
                    result.push(OptionTuple {
                        argument: Some(argument),
                        option_string: option_string.to_string(),
                        sep,
                        explicit_arg: explicit_arg.clone(),
                    });
                }
            }
        } else {
            let short_end = arg.char_indices().nth(2).map_or(arg.len(), |(at, _)| at);
            let (short_option_prefix, short_explicit_arg) = arg.split_at(short_end);
            for (argument, option_string) in self.option_strings() {
                if option_string == short_option_prefix {
                    result.push(OptionTuple {
                        argument: Some(argument),
                        option_string: option_string.to_string(),
                        sep: Some(""),
                        explicit_arg: Some(short_explicit_arg.to_string()),
                    });
                } else if option_string.starts_with(option_prefix) {
                    result.push(OptionTuple {
                        argument: Some(argument),
                        option_string: option_string.to_string(),
                        sep,
                        explicit_arg: explicit_arg.clone(),
                    });
                }
            }
        }
        result
    }

    /// The occurrences of the options on `args`, in the order they are met, or how the
    /// parser ends the run.
    pub fn parse_args(&self, args: &[OsString]) -> Result<Vec<Taken>, Exit> {
        let args: Vec<String> = args.iter().map(|arg| text_of(arg)).collect();
        Reading::new(self, &args).run()
    }
}

/// One command line while it is read.
struct Reading<'a> {
    parser: &'a ArgumentParser,
    args: &'a [String],
    /// Each argument's readings as an option, `None` for a value and for `--`.
    options: Vec<Option<Vec<OptionTuple>>>,
    /// Each argument as `O` (an option), `A` (a value) or `-` (the first `--`).
    pattern: Vec<u8>,
    /// The arguments nothing takes, in order.
    extras: Vec<String>,
    taken: Vec<Taken>,
    /// The options met so far.
    seen: BTreeSet<usize>,
}

impl<'a> Reading<'a> {
    /// Every argument is read as an option or a value before any is acted on; everything
    /// after the first `--` is a value.
    fn new(parser: &'a ArgumentParser, args: &'a [String]) -> Self {
        let mut options = Vec::with_capacity(args.len());
        let mut pattern = Vec::with_capacity(args.len());
        let mut after_dashes = false;
        for arg in args {
            if after_dashes {
                options.push(None);
                pattern.push(b'A');
            } else if arg == "--" {
                after_dashes = true;
                options.push(None);
                pattern.push(b'-');
            } else {
                let option_tuples = parser.parse_optional(arg);
                pattern.push(if option_tuples.is_some() { b'O' } else { b'A' });
                options.push(option_tuples);
            }
        }
        Reading {
            parser,
            args,
            options,
            pattern,
            extras: Vec::new(),
            taken: Vec::new(),
            seen: BTreeSet::new(),
        }
    }

    /// The parser has no positional argument: a value no option takes is left over, and
    /// the leftovers are an error only once every option was read without one.
    fn run(mut self) -> Result<Vec<Taken>, Exit> {
        let mut start_index = 0;
        if let Some(max_option_index) = self.options.iter().rposition(Option::is_some) {
            while start_index <= max_option_index {
                let next_option_index = (start_index..=max_option_index)
                    .find(|at| self.options[*at].is_some())
                    .unwrap_or(max_option_index);
                self.extras
                    .extend_from_slice(&self.args[start_index..next_option_index]);
                start_index = self.consume_optional(next_option_index)?;
            }
        }
        self.extras.extend_from_slice(&self.args[start_index..]);
        if !self.extras.is_empty() {
            return Err(Exit::Error(format!(
                "unrecognized arguments: {}",
                self.extras.join(" ")
            )));
        }
        Ok(self.taken)
    }

    /// The values the option takes of the arguments that follow it: one for an option that
    /// takes a value, which must be an argument read as a value, and none for a switch.
    fn match_argument(&self, argument: usize, pattern: &[u8]) -> Result<usize, Exit> {
        if !self.parser.arguments[argument].action.takes_a_value() {
            return Ok(0);
        }
        if pattern.first() == Some(&b'A') {
            Ok(1)
        } else {
            Err(self
                .parser
                .argument_error(argument, "expected one argument"))
        }
    }

    /// Acts on the option at `start_index` and returns the index after its value.
    fn consume_optional(&mut self, start_index: usize) -> Result<usize, Exit> {
        let option_tuples = self.options[start_index].clone().unwrap_or_default();
        if option_tuples.len() > 1 {
            let matches: Vec<&str> = option_tuples
                .iter()
                .map(|option_tuple| option_tuple.option_string.as_str())
                .collect();
            return Err(Exit::Error(format!(
                "ambiguous option: {} could match {}",
                self.args[start_index],
                matches.join(", ")
            )));
        }
        let Some(OptionTuple {
            mut argument,
            mut option_string,
            mut sep,
            mut explicit_arg,
        }) = option_tuples.into_iter().next()
        else {
            return Ok(start_index + 1);
        };
        // Several one-letter switches may stand in one argument, so one argument may act
        // on several options.
        let mut action_tuples: Vec<(usize, Vec<String>)> = Vec::new();
        let stop;
        loop {
            let Some(at) = argument else {
                self.extras.push(self.args[start_index].clone());
                return Ok(start_index + 1);
            };
            let Some(explicit) = explicit_arg.clone() else {
                let start = start_index + 1;
                let arg_count = self.match_argument(at, &self.pattern[start..])?;
                stop = start + arg_count;
                action_tuples.push((at, self.args[start..stop].to_vec()));
                break;
            };
            let arg_count = self.match_argument(at, b"A")?;
            let one_letter = !option_string[1..].starts_with('-');
            if arg_count == 0 && one_letter && !explicit.is_empty() {
                if sep == Some("=") || explicit.starts_with('-') {
                    let message = format!("ignored explicit argument {}", repr(&explicit));
                    return Err(self.parser.argument_error(at, &message));
                }
                action_tuples.push((at, Vec::new()));
                let mut chars = explicit.chars();
                option_string = format!("-{}", chars.next().unwrap_or_default());
                let rest = chars.as_str();
                match self.parser.option(&option_string) {
                    Some(next) => {
                        argument = Some(next);
                        (sep, explicit_arg) = match rest.strip_prefix('=') {
                            _ if rest.is_empty() => (None, None),
                            Some(value) => (Some("="), Some(value.to_string())),
                            None => (Some(""), Some(rest.to_string())),
                        };
                    }
                    None => {
                        self.extras.push(format!("-{explicit}"));
                        stop = start_index + 1;
                        break;
                    }
                }
            } else if arg_count == 1 {
                stop = start_index + 1;
                action_tuples.push((at, vec![explicit]));
                break;
            } else {
                let message = format!("ignored explicit argument {}", repr(&explicit));
                return Err(self.parser.argument_error(at, &message));
            }
        }
        for (at, values) in action_tuples {
            self.take_action(at, &values)?;
        }
        Ok(stop)
    }

    /// The value is converted and held to the choices before the option is held against
    /// the options it cannot stand beside.
    fn take_action(&mut self, at: usize, values: &[String]) -> Result<(), Exit> {
        let argument = &self.parser.arguments[at];
        let value = match values.first() {
            Some(text) if argument.action.takes_a_value() => self.value(at, text)?,
            _ => Value::Switch,
        };
        self.seen.insert(at);
        for group in &self.parser.exclusive {
            if !group.contains(&at) {
                continue;
            }
            if let Some(other) = group
                .iter()
                .find(|other| **other != at && self.seen.contains(other))
            {
                let message = format!(
                    "not allowed with argument {}",
                    self.parser.arguments[*other].name()
                );
                return Err(self.parser.argument_error(at, &message));
            }
        }
        if argument.action == Action::Help {
            return Err(Exit::Help);
        }
        self.taken.push(Taken {
            argument: at,
            value,
        });
        Ok(())
    }

    fn value(&self, at: usize, text: &str) -> Result<Value, Exit> {
        let argument = &self.parser.arguments[at];
        let value = match argument.kind {
            Kind::Text => Value::Text(text.to_string()),
            Kind::Path => Value::Path(path_of(text)),
            Kind::Int => match pyfmt::int_checked(text) {
                Ok(decimal) => {
                    let beyond = if decimal.starts_with('-') {
                        i64::MIN
                    } else {
                        i64::MAX
                    };
                    Value::Int(decimal.parse().unwrap_or(beyond))
                }
                Err(_) => {
                    let message = format!("invalid int value: {}", repr(text));
                    return Err(self.parser.argument_error(at, &message));
                }
            },
        };
        if let Some(choices) = argument.choices {
            if !choices.contains(&text) {
                let choices: Vec<String> = choices.iter().map(|choice| repr(choice)).collect();
                let message = format!(
                    "invalid choice: {} (choose from {})",
                    repr(text),
                    choices.join(", ")
                );
                return Err(self.parser.argument_error(at, &message));
            }
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn os(bytes: &[u8]) -> OsString {
        OsString::from_vec(bytes.to_vec())
    }

    #[test]
    fn a_byte_that_is_no_utf8_is_held_and_printed_as_the_interpreter_prints_it() {
        let text = text_of(&os(b"w\xff.txt"));
        assert_eq!(text.chars().count(), 6);
        assert_eq!(bytes_of(&text), b"w\xff.txt");
        assert_eq!(backslashreplace(&text), "w\\udcff.txt");
        assert_eq!(repr(&text), "'w\\udcff.txt'");
        // A lead byte without its continuation, a cut character and a surrogate's bytes are
        // held byte by byte.
        for bytes in [
            &b"\xe2\x28\xa1"[..],
            b"\xf0\x90\x80",
            b"\xed\xa0\x80",
            b"a\x80",
        ] {
            let text = text_of(&os(bytes));
            assert_eq!(bytes_of(&text), bytes);
            assert_eq!(text.chars().count(), bytes.len());
        }
    }

    #[test]
    fn a_character_of_the_holding_range_is_held_as_its_bytes() {
        let spelt = "a\u{10dcff}b";
        let text = text_of(OsStr::new(spelt));
        assert_eq!(text.chars().count(), 6);
        assert_eq!(bytes_of(&text), spelt.as_bytes());
        // The code point below the range is a character like any other.
        assert_eq!(text_of(OsStr::new("\u{10dc7f}")), "\u{10dc7f}");
        assert_eq!(repr("\u{10dc7f}"), "'\\U0010dc7f'");
    }

    #[test]
    fn repr_respells_the_held_bytes_alone() {
        // A backslash the text holds, before the letters of an escape, stays as written.
        let text = format!("\\U0010dcff{}", held(0xff));
        assert_eq!(repr(&text), "'\\\\U0010dcff\\udcff'");
        assert_eq!(repr("it's"), "\"it's\"");
        assert_eq!(repr("a\nb"), "'a\\nb'");
    }

    #[test]
    fn a_path_is_spelt_as_pathlib_spells_it() {
        // (the argument, `str(path)`, `path.name`, `str(path.parent)`)
        let forms = [
            ("", ".", "", "."),
            (".", ".", "", "."),
            ("./", ".", "", "."),
            ("r.json", "r.json", "r.json", "."),
            ("a/b/", "a/b", "b", "a"),
            (".//a/./b/", "a/b", "b", "a"),
            ("a/..", "a/..", "..", "a"),
            ("..", "..", "..", "."),
            ("/", "/", "", "/"),
            ("//", "//", "", "//"),
            ("///", "/", "", "/"),
            ("/a", "/a", "a", "/"),
            ("//a", "//a", "a", "//"),
            ("///a", "/a", "a", "/"),
            ("/x/../y", "/x/../y", "y", "/x/.."),
        ];
        for (argument, text, name, parent) in forms {
            let path = path_of(argument);
            assert_eq!(path_text(&path), text, "str of {argument:?}");
            assert_eq!(path_name(&path), name, "name of {argument:?}");
            assert_eq!(
                path_text(&path_parent(&path)),
                parent,
                "parent of {argument:?}"
            );
        }
        let path = path_of(&text_of(&os(b"d//w\xff.txt")));
        assert_eq!(path.as_os_str().as_bytes(), b"d/w\xff.txt");
        assert_eq!(backslashreplace(&path_name(&path)), "w\\udcff.txt");
    }

    fn parser() -> ArgumentParser {
        let mut parser = ArgumentParser::new("prog", "usage: prog\n", "help\n");
        let first = parser.add_argument("--first", Action::StoreTrue, Kind::Text, None);
        let second = parser.add_argument("--second", Action::StoreTrue, Kind::Text, None);
        parser.add_mutually_exclusive_group(&[first, second]);
        parser.add_argument("--file", Action::Store, Kind::Path, None);
        parser.add_argument("--file-list", Action::Append, Kind::Path, None);
        parser.add_argument("--kind", Action::Store, Kind::Text, Some(&["a", "b"]));
        parser.add_argument("--count", Action::Store, Kind::Int, None);
        parser
    }

    fn read(args: &[&str]) -> Result<Vec<(String, Value)>, Exit> {
        let parser = parser();
        let args: Vec<OsString> = args.iter().map(OsString::from).collect();
        Ok(parser
            .parse_args(&args)?
            .into_iter()
            .map(|taken| (parser.dest(taken.argument).to_string(), taken.value))
            .collect())
    }

    fn error(args: &[&str]) -> String {
        match read(args) {
            Err(Exit::Error(message)) => message,
            other => panic!("{args:?} is no usage error: {other:?}"),
        }
    }

    #[test]
    fn an_option_takes_its_value_from_the_next_argument_or_after_an_equals_sign() {
        let path = |text: &str| Value::Path(PathBuf::from(text));
        assert_eq!(
            read(&["--file", "a", "--file=b", "--kind=a", "--count", " +1_0 "]),
            Ok(vec![
                ("file".to_string(), path("a")),
                ("file".to_string(), path("b")),
                ("kind".to_string(), Value::Text("a".to_string())),
                ("count".to_string(), Value::Int(10)),
            ])
        );
        // An exact option string is that option though it opens a longer one.
        assert_eq!(
            read(&["--file", "x"]),
            Ok(vec![("file".to_string(), path("x"))])
        );
        assert_eq!(
            read(&["--file-l", "x", "--fi=y=z"]).unwrap_err(),
            Exit::Error(
                "ambiguous option: --fi=y=z could match --first, --file, --file-list".to_string()
            )
        );
        assert_eq!(
            read(&["--k", "b"]),
            Ok(vec![("kind".to_string(), Value::Text("b".to_string()))])
        );
        // A value that opens with a dash is taken where it reads as a negative number or
        // holds a space, and after an equals sign.
        assert_eq!(
            read(&["--count", "-1"]),
            Ok(vec![("count".to_string(), Value::Int(-1))])
        );
        assert_eq!(
            read(&["--file", "-x y"]),
            Ok(vec![("file".to_string(), path("-x y"))])
        );
        assert_eq!(
            read(&["--file=-h"]),
            Ok(vec![("file".to_string(), path("-h"))])
        );
        assert_eq!(
            read(&["--file", "-"]),
            Ok(vec![("file".to_string(), path("-"))])
        );
        assert_eq!(
            read(&["--count", "99999999999999999999"]),
            Ok(vec![("count".to_string(), Value::Int(i64::MAX))])
        );
    }

    #[test]
    fn usage_errors_carry_the_interpreters_texts() {
        assert_eq!(error(&["--file"]), "argument --file: expected one argument");
        assert_eq!(
            error(&["--file", "-x"]),
            "argument --file: expected one argument"
        );
        assert_eq!(
            error(&["--file", "--"]),
            "argument --file: expected one argument"
        );
        assert_eq!(
            error(&["--file", "--help"]),
            "argument --file: expected one argument"
        );
        assert_eq!(
            error(&["--kind", "c"]),
            "argument --kind: invalid choice: 'c' (choose from 'a', 'b')"
        );
        assert_eq!(
            error(&["--count", "x"]),
            "argument --count: invalid int value: 'x'"
        );
        assert_eq!(
            error(&["--count=1.5"]),
            "argument --count: invalid int value: '1.5'"
        );
        assert_eq!(
            error(&["--first", "--second"]),
            "argument --second: not allowed with argument --first"
        );
        assert_eq!(
            error(&["--first=1"]),
            "argument --first: ignored explicit argument '1'"
        );
        assert_eq!(
            error(&["--help="]),
            "argument -h/--help: ignored explicit argument ''"
        );
        assert_eq!(
            error(&["-h=1"]),
            "argument -h/--help: ignored explicit argument '1'"
        );
        assert_eq!(
            error(&["--f", "x"]),
            "ambiguous option: --f could match --first, --file, --file-list"
        );
        assert_eq!(
            error(&["stray", "--first", "--", "--file"]),
            "unrecognized arguments: stray -- --file"
        );
        assert_eq!(
            error(&["--bogus", "value", "--first"]),
            "unrecognized arguments: --bogus value"
        );
        // The first error ends the reading; an unknown option is an error only after every
        // other argument was read.
        assert_eq!(
            error(&["--bogus", "--kind", "c", "--help"]),
            "argument --kind: invalid choice: 'c' (choose from 'a', 'b')"
        );
    }

    #[test]
    fn help_ends_the_reading_where_it_is_met() {
        for args in [
            &["-h"][..],
            &["--help"],
            &["--he"],
            &["-hx"],
            &["-help"],
            &["-hh"],
        ] {
            assert_eq!(read(args), Err(Exit::Help), "{args:?}");
        }
        assert_eq!(read(&["--bogus", "--help"]), Err(Exit::Help));
        assert_eq!(read(&["--help", "--kind", "c"]), Err(Exit::Help));
        assert_eq!(
            parser().error_text("unrecognized arguments: x"),
            "usage: prog\nprog: error: unrecognized arguments: x\n"
        );
    }
}
