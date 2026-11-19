// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Perl's regular expression dialect, as `/$pattern/i` reads a filter's
//! pattern in Ensembl VEP (`FilterSet.pm` lines 714-715), rewritten for the
//! `regex` crate.
//!
//! [`translate`] reads the bytes of a pattern twice, a first pass for the
//! rule set and the group count and a second that writes the pattern the
//! crate compiles, over one of two readings ([`Rules`]):
//!
//! * Under Perl's default rules over a byte string a byte above ASCII has no
//!   case and is no `\w`, `\s` or `\d`, so the pattern runs over the bytes of
//!   the value with the crate's Unicode mode off.
//! * A `\p{}` or `\P{}` class, a `\x{}`, `\o{}` or octal escape above 255 or
//!   a `\N{U+}` escape anywhere switches the whole pattern to Perl's Unicode
//!   rules, under which each byte is the Latin-1 character of its value; the
//!   pattern then runs over that reading with the crate's Unicode mode on. A
//!   `u`, `a` or `aa` charset flag switches its own group, or the rest of its
//!   group, alone: the pattern runs over the Latin-1 reading, and where the
//!   default `d` charset stays in force the byte rules are written out, with
//!   ASCII classes, no case above ASCII and no fold to a sequence.
//!
//! Under Unicode rules a case-insensitive match compares full case folds:
//! `ß` and `ẞ` fold to `ss`, the ligatures `ﬀ ﬁ ﬂ ﬃ ﬄ ﬅ ﬆ` to `ff fi fl ffi
//! ffl st st`, and the fold of a run of literal characters is compared whole,
//! so `ss` matches a value's `ß` and `\x{DF}` matches `ss`, `SS` and `sS`. A
//! run is a maximal sequence of adjacent literals in one branch under one
//! charset. An unquantified non-capturing group of one branch or of only
//! empty branches, a flag group that keeps case sensitivity, a comment and a
//! bracketed class whose members all fold to one character join it; a
//! quantifier, a capturing group, a group of several branches unless all are
//! empty, a class of several fold sets, an assertion and a change of case
//! sensitivity or of charset break it. A group that sets or clears `m` or
//! `s`, by its opening or by a flag group inside, joins as well: the join
//! drops those flags and writes them back around each atom inside that it
//! writes as text, since only `.`, `^` and `$` read them. A quantified literal
//! with a full fold repeats the character or its fold string; a bracketed
//! class matches the fold string of each listed member with such a fold,
//! never of a range, and a negated class never does. Under `aa` no fold
//! crosses ASCII, so `ß` stays `ß`; under `a` and `aa` `\w`, `\d`, `\s`, `\b`
//! and the POSIX classes are ASCII. Case-insensitive, `[[:lower:]]` and
//! `[[:upper:]]` are the cased letters.
//!
//! The brace and escape forms follow Perl: `{` opens a quantifier only as
//! `{n}`, `{n,}`, `{n,m}` or `{,m}`, blanks allowed inside, after an atom, and
//! is a literal elsewhere; a leading zero or a bound above 65534 is refused;
//! `{n,m}` with `n` above `m` never matches. `\o{}` and the octal escapes `\0`
//! to `\777` are characters by value, `\1` to `\9` are back-references, and
//! `\10` and above is a back-reference when that many groups open before it
//! and octal otherwise. `\N{U+}` holds hexadecimal numbers, single
//! underscores allowed between the digits, separated by `.` into a sequence of
//! characters that is one atom, a run of literals of its own and an
//! alternative inside a bracketed class. `\e`, `\a`, `\cX`, `\h`, `\v`, `\R`,
//! `\N`, `\K`, `\G` and `\Z` are read as Perl reads them, and an escape Perl
//! does not know is its character. Back-references, look-around and the other
//! constructs the crate lacks are refused naming the pattern.

use std::fmt::Write as _;

use crate::filter_input;

/// The reading a pattern runs over.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) enum Rules {
    /// Perl's default rules over a byte string: the crate's Unicode mode off.
    #[default]
    Bytes,
    /// The bytes read as Latin-1 characters with the crate's Unicode mode
    /// on: Perl's Unicode rules, and its byte rules where a `d` charset is in
    /// force in a pattern without a Unicode trigger.
    Latin1,
}

/// The rules of the whole pattern, against which a point's charset flags
/// settle the rules in force there.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct RuleSet {
    rules: Rules,
    /// Whether a Unicode trigger puts every `d` scope under Unicode rules.
    wide_unicode: bool,
}

impl RuleSet {
    /// Whether Perl's byte rules hold at a point of the Latin-1 reading: a
    /// `d` charset in a pattern with a charset flag but no Unicode trigger.
    /// Over the Latin-1 reading they are written out: ASCII classes and no
    /// case above ASCII.
    fn bytes_over_latin1(self, flags: Flags) -> bool {
        self.rules == Rules::Latin1 && !self.wide_unicode && !flags.unicode && flags.ascii == 0
    }

    /// Whether `\w`, `\d`, `\s`, `\b` and the POSIX classes are written as
    /// ASCII classes: under `a` and `aa`, and under byte rules over the
    /// Latin-1 reading. Over the bytes the crate's own classes are ASCII.
    fn ascii_classes(self, flags: Flags) -> bool {
        self.rules == Rules::Latin1 && (flags.ascii > 0 || self.bytes_over_latin1(flags))
    }

    /// Whether literals fold fully: Unicode rules in force, case-insensitive,
    /// no `aa`.
    fn folds(self, flags: Flags) -> bool {
        self.rules == Rules::Latin1 && flags.i && flags.ascii != 2 && !self.bytes_over_latin1(flags)
    }

    /// How a run of literals is written under `flags`: case sensitivity,
    /// whether it folds fully, and whether `aa` keeps its characters that
    /// fold to ASCII case-sensitive. A literal under another state ends the
    /// run, so `u` and `a` scopes, both folding, share one and a `d` scope
    /// joins them only where a Unicode trigger puts it under Unicode rules.
    fn fold_state(self, flags: Flags) -> (bool, bool, bool) {
        (flags.i, self.folds(flags), flags.ascii == 2)
    }
}

/// Why a pattern does not compile.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// Perl's own refusal, with its message (`Invalid quantifier in {,}`).
    Perl(String),
    /// A construct the supported dialect lacks.
    Unsupported(String),
}

/// A pattern rewritten for the crate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Translated {
    pub(crate) rules: Rules,
    pub(crate) pattern: String,
}

/// `pattern` as the crate reads it, under the rules Perl would read it.
pub(crate) fn translate(pattern: &[u8]) -> Result<Translated, Refusal> {
    let flags = Flags::default();
    // The first reading finds the rule set and the group count, which the
    // digit escapes of the second depend on.
    // A refusal of the first reading recurs in the second, which reports it.
    let mut survey = Reader::new(pattern, RuleSet::default(), None);
    let _ = survey.branches(flags);
    // A Unicode trigger puts the whole pattern under Unicode rules; a
    // charset flag alone puts its scope there and leaves the rest under
    // byte rules, written over the same Latin-1 reading.
    let rule_set = RuleSet {
        rules: if survey.unicode || survey.charset_flag {
            Rules::Latin1
        } else {
            Rules::Bytes
        },
        wide_unicode: survey.unicode,
    };
    let mut reader = Reader::new(pattern, rule_set, Some(survey.groups));
    let branches = match reader.branches(flags) {
        Ok(branches) if reader.at == pattern.len() => branches,
        Ok(_) => return Err(Refusal::Perl("Unmatched )".to_string())),
        // Perl reads past a construct the crate lacks, so its own refusal
        // of a missing group comes first. The first reading stopped at the
        // construct too, so the groups it never counted are at most the `(`
        // bytes from there on.
        Err(Refusal::Unsupported(message)) => {
            let uncounted = pattern[survey.at..].iter().filter(|b| **b == b'(').count();
            return Err(match reader.missing_reference {
                Some(number) if number > survey.groups + uncounted => perl(MISSING_GROUP),
                _ => Refusal::Unsupported(message),
            });
        }
        Err(refusal) => return Err(refusal),
    };
    if reader.missing_reference.is_some() {
        return Err(perl(MISSING_GROUP));
    }
    let mut emitter = Emitter {
        rule_set,
        out: String::with_capacity(pattern.len()),
        longest_run: 0,
    };
    emitter.branches(&branches, flags);
    if emitter.longest_run > LONGEST_RUN {
        return Err(Refusal::Unsupported(format!(
            "a run of more than {LONGEST_RUN} letters s under full case folding"
        )));
    }
    Ok(Translated {
        rules: rule_set.rules,
        pattern: emitter.out,
    })
}

const LOOK_AROUND: &str = "look-around, including look-ahead and look-behind, is not supported";

const MISSING_GROUP: &str = "Reference to nonexistent group";

/// The flags in force at a point of the pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Flags {
    /// Case-insensitive: on at the start, since the pattern runs under `/i`.
    i: bool,
    /// Blanks and comments ignored: 0 off, 1 `x`, 2 `xx`.
    x: u8,
    /// ASCII restriction: 0 off, 1 `a`, 2 `aa`.
    ascii: u8,
    /// `u`: Unicode rules in this scope whatever the rest of the pattern
    /// reads under; off with `ascii` at 0 is Perl's default `d` charset.
    unicode: bool,
    /// `n`: plain groups do not capture.
    n: bool,
}

impl Default for Flags {
    fn default() -> Flags {
        Flags {
            i: true,
            x: 0,
            ascii: 0,
            unicode: false,
            n: false,
        }
    }
}

#[derive(Debug)]
enum Tok {
    /// A literal character.
    Lit(char),
    Class(Class),
    Group(Group),
    /// A flag group `(?i)`, with the text the crate reads (empty when the
    /// crate needs none) and the flags after it.
    Flags {
        text: String,
        flags: Flags,
    },
    /// The quantifier of the token before it; `never` for `{n,m}` with `n`
    /// above `m`, which never matches.
    Quant {
        text: String,
        never: bool,
    },
    /// Text the crate reads as written, one atom.
    Other(String),
}

#[derive(Debug)]
struct Class {
    /// The bracketed class as the crate reads it; empty when the class lists
    /// `\N{U+}` sequences alone.
    text: String,
    negated: bool,
    /// Every single character listed, not from a range.
    members: Vec<char>,
    /// The `\N{U+}` sequences of several characters listed, each an
    /// alternative to the class.
    sequences: Vec<Vec<char>>,
    /// Whether every item is a single character of one fold set.
    uniform: bool,
}

#[derive(Debug)]
struct Group {
    capturing: bool,
    /// The opening text the crate reads: `(`, `(?P<name>`, `(?:`, `(?i:`.
    open: String,
    branches: Vec<Vec<Tok>>,
    /// The flags inside.
    flags: Flags,
}

/// One item of a bracketed class.
enum Item {
    Char(char),
    /// A class escape or a POSIX class, as the crate reads it.
    Set(String),
    /// A `\N{U+}` escape of several characters.
    Sequence(Vec<char>),
}

struct Reader<'p> {
    pattern: &'p [u8],
    at: usize,
    rule_set: RuleSet,
    /// The capturing groups of the whole pattern; `None` on the first
    /// reading, which counts them.
    captures: Option<usize>,
    /// Whether a Unicode trigger was read, which puts the whole pattern
    /// under Unicode rules: a `\p{}` or `\P{}` class, a `\x{}`, `\o{}` or
    /// octal escape above 255 or a `\N{U+}` escape.
    unicode: bool,
    /// Whether a `u`, `a` or `aa` charset flag was read, which puts its own
    /// scope under Unicode rules.
    charset_flag: bool,
    /// The capturing groups opened so far.
    groups: usize,
    /// The open groups at this point.
    depth: usize,
    /// The largest number of a reference to a group beyond `captures`, which
    /// Perl refuses once the whole pattern is read.
    missing_reference: Option<usize>,
}

/// The deepest group nesting read: a pattern read in a thread of the
/// default stack size stays well inside it.
const NEST_LIMIT: usize = 100;

/// The blanks `x` skips: Perl's `\s` on a byte string.
fn is_blank(byte: u8) -> bool {
    filter_input::perl_space(char::from(byte))
}

fn perl(message: &str) -> Refusal {
    Refusal::Perl(message.to_string())
}

impl<'p> Reader<'p> {
    fn new(pattern: &'p [u8], rule_set: RuleSet, captures: Option<usize>) -> Reader<'p> {
        Reader {
            pattern,
            at: 0,
            rule_set,
            captures,
            unicode: false,
            charset_flag: false,
            groups: 0,
            depth: 0,
            missing_reference: None,
        }
    }

    fn peek(&self) -> Option<u8> {
        self.pattern.get(self.at).copied()
    }

    fn peek_at(&self, offset: usize) -> Option<u8> {
        self.pattern.get(self.at + offset).copied()
    }

    /// Under `x`, past blanks and `#` comments.
    fn skip_blanks(&mut self, flags: Flags) {
        if flags.x == 0 {
            return;
        }
        while let Some(byte) = self.peek() {
            if is_blank(byte) {
                self.at += 1;
            } else if byte == b'#' {
                while self.peek().is_some_and(|b| b != b'\n') {
                    self.at += 1;
                }
            } else {
                break;
            }
        }
    }

    /// Past what Perl skips between an atom, its quantifier and the
    /// quantifier's modifier: the blanks and `#` comments of `x`, and a
    /// `(?#...)` comment under any flags.
    fn skip_ignored(&mut self, flags: Flags) -> Result<(), Refusal> {
        loop {
            self.skip_blanks(flags);
            if !self.pattern[self.at..].starts_with(b"(?#") {
                return Ok(());
            }
            self.comment()?;
        }
    }

    /// Past the `(?#...)` comment that `self.at` lies in, which ends at the
    /// first `)`.
    fn comment(&mut self) -> Result<(), Refusal> {
        let Some(close) = self.pattern[self.at..].iter().position(|b| *b == b')') else {
            return Err(perl("Sequence (?#... not terminated"));
        };
        self.at += close + 1;
        Ok(())
    }

    /// The branches up to a `)`, left unread, or the end.
    fn branches(&mut self, mut flags: Flags) -> Result<Vec<Vec<Tok>>, Refusal> {
        let mut branches = vec![Vec::new()];
        loop {
            self.skip_blanks(flags);
            let Some(byte) = self.peek() else { break };
            match byte {
                b')' => break,
                b'|' => {
                    self.at += 1;
                    branches.push(Vec::new());
                }
                _ => {
                    let branch = branches.last_mut().expect("one branch");
                    self.item(&mut flags, branch)?;
                }
            }
        }
        Ok(branches)
    }

    /// One atom with its quantifier, a flag group or a comment.
    fn item(&mut self, flags: &mut Flags, branch: &mut Vec<Tok>) -> Result<(), Refusal> {
        let byte = self.peek().expect("a byte");
        let tok = match byte {
            b'\\' => {
                self.at += 1;
                self.escape(*flags)?
            }
            b'[' => {
                self.at += 1;
                self.class(*flags)?
            }
            b'(' => {
                self.at += 1;
                match self.group(flags)? {
                    Some(tok) => tok,
                    None => return Ok(()),
                }
            }
            // A `{` that follows no atom is a literal, a quantifier shape or
            // not; after an atom the quantifier was read with the atom, so
            // this `{` opens none.
            b'{' | b'}' => {
                self.at += 1;
                Tok::Lit(char::from(byte))
            }
            b'*' | b'+' | b'?' => return Err(perl("Quantifier follows nothing")),
            b'.' | b'^' | b'$' => {
                self.at += 1;
                Tok::Other(char::from(byte).to_string())
            }
            _ => {
                self.at += 1;
                Tok::Lit(char::from(byte))
            }
        };
        let is_atom = !matches!(tok, Tok::Flags { .. });
        branch.push(tok);
        if is_atom {
            if let Some(quant) = self.quantifier(*flags)? {
                branch.push(quant);
            }
        }
        Ok(())
    }

    /// The quantifier after an atom, with a lazy `?` kept and a possessive
    /// `+` read as greedy; `None` when none follows. A quantifier on a
    /// quantifier is Perl's `Nested quantifiers` refusal.
    fn quantifier(&mut self, flags: Flags) -> Result<Option<Tok>, Refusal> {
        self.skip_ignored(flags)?;
        let (mut text, never) = match self.peek() {
            Some(byte @ (b'*' | b'+' | b'?')) => {
                self.at += 1;
                (char::from(byte).to_string(), false)
            }
            Some(b'{') => match self.brace(self.at) {
                Some(Ok((text, never, end))) => {
                    self.at = end;
                    (text, never)
                }
                Some(Err(refusal)) => return Err(refusal),
                None => return Ok(None),
            },
            _ => return Ok(None),
        };
        self.skip_ignored(flags)?;
        match self.peek() {
            Some(b'?') => {
                self.at += 1;
                text.push('?');
            }
            Some(b'+') => self.at += 1,
            _ => {}
        }
        self.skip_ignored(flags)?;
        let nested = match self.peek() {
            Some(b'*' | b'+' | b'?') => true,
            Some(b'{') => self.brace(self.at).is_some(),
            _ => false,
        };
        if nested {
            return Err(perl("Nested quantifiers"));
        }
        Ok(Some(Tok::Quant { text, never }))
    }

    /// The brace quantifier whose `{` is `pattern[at]`: `{n}`, `{n,}`,
    /// `{n,m}` or `{,m}` with blanks allowed inside, as the crate's text, its
    /// never-matching flag and the index after the `}`; `None` when the brace
    /// opens no quantifier; Perl's refusal for a leading zero or a bound
    /// above 65534.
    fn brace(&self, at: usize) -> Option<Result<(String, bool, usize), Refusal>> {
        let p = self.pattern;
        let blanks = |mut i: usize| {
            while p.get(i).is_some_and(|b| matches!(b, b' ' | b'\t')) {
                i += 1;
            }
            i
        };
        let digits = |i: usize| {
            let end = i + p[i..].iter().take_while(|b| b.is_ascii_digit()).count();
            (i < end).then(|| &p[i..end])
        };
        let bound = |d: &[u8]| -> Result<u32, Refusal> {
            if d.len() > 1 && d[0] == b'0' {
                return Err(perl("Invalid quantifier in {,}"));
            }
            let value = d.iter().fold(0u32, |v, b| {
                v.saturating_mul(10).saturating_add((b - b'0').into())
            });
            if value > 65534 {
                return Err(perl("Quantifier in {,} bigger than 65534"));
            }
            Ok(value)
        };
        let mut i = blanks(at + 1);
        let low = digits(i);
        i = blanks(i + low.map_or(0, <[u8]>::len));
        match p.get(i) {
            Some(b'}') => {
                let low = low?;
                Some(bound(low).map(|n| (format!("{{{n}}}"), false, i + 1)))
            }
            Some(b',') => {
                i = blanks(i + 1);
                let high = digits(i);
                i = blanks(i + high.map_or(0, <[u8]>::len));
                if p.get(i) != Some(&b'}') || (low.is_none() && high.is_none()) {
                    return None;
                }
                let range = || -> Result<(String, bool, usize), Refusal> {
                    let low = low.map(bound).transpose()?.unwrap_or(0);
                    let high = high.map(bound).transpose()?;
                    let text = match high {
                        Some(high) => format!("{{{low},{high}}}"),
                        None => format!("{{{low},}}"),
                    };
                    Ok((text, high.is_some_and(|h| h < low), i + 1))
                };
                Some(range())
            }
            _ => None,
        }
    }

    /// The escape after a backslash.
    fn escape(&mut self, flags: Flags) -> Result<Tok, Refusal> {
        let Some(c) = self.peek() else {
            return Err(perl("Trailing \\"));
        };
        self.at += 1;
        Ok(match c {
            b'x' => Tok::Lit(self.hex()?),
            b'o' => Tok::Lit(self.braced_octal()?),
            b'0' => {
                self.at -= 1;
                Tok::Lit(self.octal())
            }
            b'1'..=b'9' => self.digits(c)?,
            b'N' => self.named(flags)?,
            b'p' | b'P' => Tok::Other(self.property(c)?),
            b'e' => Tok::Lit('\x1B'),
            b'a' => Tok::Lit('\x07'),
            b't' => Tok::Lit('\t'),
            b'n' => Tok::Lit('\n'),
            b'r' => Tok::Lit('\r'),
            b'f' => Tok::Lit('\x0C'),
            b'c' => Tok::Lit(self.control()?),
            b'd' | b'D' | b'w' | b'W' | b's' | b'S' | b'h' | b'H' | b'v' | b'V' => {
                Tok::Other(class_escape(c, flags, self.rule_set))
            }
            b'R' => Tok::Other(r"(?:\r\n|[\n\x0B\x0C\r\x85])".to_string()),
            b'b' | b'B' if self.peek() == Some(b'{') => return Err(self.bound_type(c)),
            b'b' | b'B' => Tok::Other(word_boundary(c, flags, self.rule_set)),
            b'A' | b'G' => Tok::Other(r"\A".to_string()),
            b'z' | b'Z' => Tok::Other(r"\z".to_string()),
            b'K' => Tok::Other("(?:)".to_string()),
            // The crate refuses these: back-references, `\X`, `\C`.
            b'k' | b'g' | b'X' | b'C' => Tok::Other(format!("\\{}", char::from(c))),
            // An escape Perl does not know is its character; a backslash
            // before a byte above ASCII is that byte.
            _ => Tok::Lit(char::from(c)),
        })
    }

    /// `\x`: the braced or unbraced hex escape after the `x`.
    fn hex(&mut self) -> Result<char, Refusal> {
        if self.peek() == Some(b'{') {
            let Some(braced) = braced_number(self.pattern, self.at, 16) else {
                return Err(perl("Missing right brace on \\x{}"));
            };
            return self.braced(braced);
        }
        let (code, after) = unbraced_hex(self.pattern, self.at);
        self.at = after;
        Ok(scalar(code))
    }

    /// The character of a braced escape: a value above Perl's largest code
    /// point, 0x7FFFFFFFFFFFFFFF, is Perl's refusal, which names the value in
    /// hex below 0xFFFFFFFFFFFFFFFF, Perl's overflow mark, and the digits
    /// read in the escape's own radix from there; one above 255 switches the
    /// rules.
    fn braced(&mut self, braced: Braced) -> Result<char, Refusal> {
        let digits = &self.pattern[braced.digits.clone()];
        let code = code_point(braced.value, digits, braced.radix, false)?;
        self.at = braced.after;
        Ok(self.character(code))
    }

    /// The character of a code read from an escape, recording a code above
    /// 255 as a switch to Unicode rules.
    fn character(&mut self, code: u32) -> char {
        if code > 0xFF {
            self.unicode = true;
        }
        scalar(code)
    }

    /// `\o{}`: the braced octal escape after the `o`.
    fn braced_octal(&mut self) -> Result<char, Refusal> {
        if self.peek() != Some(b'{') {
            return Err(perl("Missing braces on \\o{}"));
        }
        let Some(braced) = braced_number(self.pattern, self.at, 8) else {
            return Err(perl("Missing right brace on \\o{}"));
        };
        if braced.empty {
            return Err(perl("Empty \\o{}"));
        }
        self.braced(braced)
    }

    /// Up to three octal digits at `self.at`, as Perl reads `\0` and the
    /// octal reading of `\1` to `\7`.
    fn octal(&mut self) -> char {
        let digits = self.pattern[self.at..]
            .iter()
            .take(3)
            .take_while(|b| (b'0'..=b'7').contains(b))
            .count();
        let value = self.pattern[self.at..self.at + digits]
            .iter()
            .fold(0u32, |v, b| v * 8 + u32::from(b - b'0'));
        self.at += digits;
        self.character(value)
    }

    /// `\1` to `\9` and the longer digit runs: `\1` to `\9` are
    /// back-references, as is a longer number not above the count of groups
    /// opened before it; a longer number above that count is octal when it
    /// starts with `1` to `7`, and a reference above the whole pattern's
    /// count is Perl's refusal.
    fn digits(&mut self, first: u8) -> Result<Tok, Refusal> {
        let start = self.at - 1;
        let digits = self.pattern[start..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .count();
        let number = self.pattern[start..start + digits]
            .iter()
            .fold(0usize, |v, b| {
                v.saturating_mul(10).saturating_add((b - b'0').into())
            });
        if digits > 1 && first < b'8' && number > self.groups {
            self.at = start;
            return Ok(Tok::Lit(self.octal()));
        }
        // The first reading has no group count yet and judges no reference.
        if self.captures.is_some_and(|captures| number > captures) {
            self.missing_reference = self.missing_reference.max(Some(number));
        }
        self.at = start + digits;
        Ok(Tok::Other(format!("\\{number}")))
    }

    /// `\N{U+hex}` is the character and `\N{U+hex.hex}` the sequence of
    /// characters, one atom, read by [`Reader::code_points`]; `\N{name}` is
    /// refused; `\N` is a non-newline.
    fn named(&mut self, flags: Flags) -> Result<Tok, Refusal> {
        if self.peek() != Some(b'{') {
            return Ok(Tok::Other(r"[^\n]".to_string()));
        }
        let mut chars = self.code_points()?;
        Ok(if chars.len() == 1 {
            Tok::Lit(chars.remove(0))
        } else {
            Tok::Group(Group {
                capturing: false,
                open: "(?:".to_string(),
                branches: vec![chars.into_iter().map(Tok::Lit).collect()],
                flags,
            })
        })
    }

    /// The characters of the `\N{U+...}` escape whose `{` is at `self.at`:
    /// hexadecimal numbers, single underscores allowed between the digits of
    /// one, separated by `.`; anything else is Perl's refusal, as is an empty
    /// number. Every `\N{U+}` escape switches the pattern to Unicode rules.
    fn code_points(&mut self) -> Result<Vec<char>, Refusal> {
        let Some(close) = self.pattern[self.at..].iter().position(|b| *b == b'}') else {
            return Err(perl("Missing right brace on \\N{}"));
        };
        let content = &self.pattern[self.at + 1..self.at + close];
        self.at += close + 1;
        self.unicode = true;
        let Some(numbers) = content.strip_prefix(b"U+") else {
            return Err(Refusal::Unsupported(format!(
                "the escape \\N{{{}}} names a character",
                String::from_utf8_lossy(content)
            )));
        };
        let invalid = || perl("Invalid hexadecimal number in \\N{U+...}");
        numbers
            .split(|b| *b == b'.')
            .map(|digits| {
                let well_formed = !digits.is_empty()
                    && digits.first().is_some_and(u8::is_ascii_hexdigit)
                    && digits.last().is_some_and(u8::is_ascii_hexdigit)
                    && digits
                        .windows(2)
                        .all(|w| w[0].is_ascii_hexdigit() || w[1].is_ascii_hexdigit())
                    && digits.iter().all(|b| b.is_ascii_hexdigit() || *b == b'_');
                if !well_formed {
                    return Err(invalid());
                }
                let value = digits
                    .iter()
                    .filter_map(|b| char::from(*b).to_digit(16))
                    .fold(0u128, |v, d| {
                        v.saturating_mul(16).saturating_add(u128::from(d))
                    });
                code_point(value, digits, 16, true).map(scalar)
            })
            .collect()
    }

    /// `\p{..}`, `\P{..}`, `\pL`, as written.
    fn property(&mut self, c: u8) -> Result<String, Refusal> {
        self.unicode = true;
        let start = self.at;
        if self.peek() == Some(b'{') {
            let Some(close) = self.pattern[start..].iter().position(|b| *b == b'}') else {
                return Err(perl("Missing right brace on \\p{}"));
            };
            self.at = start + close + 1;
        } else if self.peek().is_some() {
            self.at += 1;
        } else {
            return Err(perl("Empty \\p"));
        }
        Ok(format!(
            "\\{}{}",
            char::from(c),
            String::from_utf8_lossy(&self.pattern[start..self.at])
        ))
    }

    /// `\b{...}` or `\B{...}` after the letter: a bound type, which the crate
    /// lacks, or a name Perl refuses, read as written.
    fn bound_type(&self, c: u8) -> Refusal {
        let letter = char::from(c);
        let Some(close) = self.pattern[self.at..].iter().position(|b| *b == b'}') else {
            return Refusal::Perl(format!("Missing right brace on \\{letter}{{}}"));
        };
        let name = String::from_utf8_lossy(&self.pattern[self.at + 1..self.at + close]);
        if name.is_empty() {
            return Refusal::Perl(format!("Empty \\{letter}{{}}"));
        }
        if matches!(name.as_ref(), "wb" | "gcb" | "lb" | "sb" | "g") {
            return Refusal::Unsupported(format!(
                "the bound type \\{letter}{{{name}}} is not supported"
            ));
        }
        Refusal::Perl(format!("'{name}' is an unknown bound type"))
    }

    /// `\cX`: the control character of `X`, which Perl requires to be
    /// printable ASCII other than `{`.
    fn control(&mut self) -> Result<char, Refusal> {
        let Some(x) = self.peek().filter(|x| (0x20..=0x7E).contains(x)) else {
            return Err(perl("Character following \"\\c\" must be printable ASCII"));
        };
        if x == b'{' {
            return Err(perl("Use \";\" instead of \"\\c{\""));
        }
        self.at += 1;
        Ok(char::from(x.to_ascii_uppercase() ^ 0x40))
    }

    /// The index past the spaces and tabs at `at`, which `xx` ignores inside
    /// a bracketed class; `at` itself under other flags.
    fn past_class_blanks(&self, at: usize, flags: Flags) -> usize {
        if flags.x != 2 {
            return at;
        }
        at + self.pattern[at..]
            .iter()
            .take_while(|b| matches!(b, b' ' | b'\t'))
            .count()
    }

    /// Where the end of a range starts when the `-` at `dash` opens one: past
    /// the `-` and the blanks `xx` skips. `None` when no `-` is at `dash` or a
    /// `]` follows it, which makes the `-` a member.
    fn range_end(&self, dash: usize, flags: Flags) -> Option<usize> {
        if self.pattern.get(dash) != Some(&b'-') {
            return None;
        }
        let end = self.past_class_blanks(dash + 1, flags);
        let opens_range = self.pattern.get(end).is_some_and(|b| *b != b']');
        opens_range.then_some(end)
    }

    /// A bracketed class after its `[`. Under `xx` a space or tab is skipped
    /// before the `^`, before the leading `]` that is a member, before each
    /// item and between a range's `-` and its end.
    fn class(&mut self, flags: Flags) -> Result<Tok, Refusal> {
        self.at = self.past_class_blanks(self.at, flags);
        let negated = self.peek() == Some(b'^');
        if negated {
            self.at += 1;
        }
        let mut items: Vec<ClassItem> = Vec::new();
        let mut sequences = Vec::new();
        // The dialect's refusals of a sequence wait for the `]`, so that Perl's
        // own refusal of an unterminated class comes first.
        let mut sequence_ends_range = false;
        let mut sequence_before_dash = false;
        let mut first = true;
        // The last single character read, which a `-` after it ranges from.
        let mut pending: Option<char> = None;
        loop {
            self.at = self.past_class_blanks(self.at, flags);
            let Some(byte) = self.peek() else {
                return Err(perl("Unmatched ["));
            };
            if byte == b']' && !first {
                self.at += 1;
                break;
            }
            first = false;
            // A `-` after a character opens a range unless `]` follows it.
            let range = if pending.is_some() {
                self.range_end(self.at, flags)
            } else {
                None
            };
            if let Some(end_at) = range {
                self.at = end_at;
                let start = pending.take().expect("a range start");
                match self.class_item(flags)? {
                    Item::Char(end) => {
                        if start > end {
                            return Err(perl("Invalid [] range"));
                        }
                        if start != end {
                            items.pop();
                            items.push(ClassItem::Range(start, end));
                        }
                    }
                    Item::Set(set) => {
                        items.push(ClassItem::Char('-'));
                        items.push(ClassItem::Set(set));
                    }
                    Item::Sequence(_) => sequence_ends_range = true,
                }
                continue;
            }
            match self.class_item(flags)? {
                Item::Char(ch) => {
                    items.push(ClassItem::Char(ch));
                    pending = Some(ch);
                }
                Item::Set(set) => {
                    items.push(ClassItem::Set(set));
                    pending = None;
                }
                Item::Sequence(chars) => {
                    // Perl reads a `-` after a sequence oddly, as a member or
                    // as a range from the letter `N` of the escape to the next
                    // item, refused as inverted below `N`; the item is read
                    // here and dropped so that the rest of the class is read
                    // as Perl reads it.
                    let dash = self.past_class_blanks(self.at, flags);
                    if self.pattern.get(dash) == Some(&b'-') {
                        sequence_before_dash = true;
                        if let Some(end_at) = self.range_end(dash, flags) {
                            self.at = end_at;
                            if let Item::Char(end) = self.class_item(flags)? {
                                if end < 'N' {
                                    return Err(perl("Invalid [] range"));
                                }
                            }
                        }
                    }
                    sequences.push(chars);
                    pending = None;
                }
            }
        }
        if sequence_ends_range {
            return Err(sequence_in_class("as the end of a range"));
        }
        if sequence_before_dash {
            return Err(sequence_in_class("before a `-`"));
        }
        if negated && !sequences.is_empty() {
            return Err(sequence_in_class("in a negated class"));
        }
        let members: Vec<char> = items
            .iter()
            .filter_map(|item| match item {
                ClassItem::Char(ch) => Some(*ch),
                _ => None,
            })
            .collect();
        let uniform = !negated
            && sequences.is_empty()
            && items.len() == members.len()
            && members.first().is_some_and(|&first| {
                let key = fold_key(first);
                members.iter().all(|&member| fold_key(member) == key)
            });
        let text = if items.is_empty() {
            String::new()
        } else {
            class_text(&items, negated, flags, self.rule_set)
        };
        Ok(Tok::Class(Class {
            text,
            negated,
            members,
            sequences,
            uniform,
        }))
    }

    /// One item of a bracketed class at `self.at`.
    fn class_item(&mut self, flags: Flags) -> Result<Item, Refusal> {
        let byte = self.peek().expect("a byte");
        self.at += 1;
        if byte == b'[' {
            if self.peek() == Some(b':') {
                let rest = &self.pattern[self.at + 1..];
                if let Some(end) = rest.windows(2).position(|w| w == b":]") {
                    let span = &rest[..end];
                    let written = String::from_utf8_lossy(span).into_owned();
                    let (negated, name) = match written.strip_prefix('^') {
                        Some(name) => (true, name),
                        None => (false, written.as_str()),
                    };
                    if let Some(set) = posix_class(name, negated, flags, self.rule_set) {
                        self.at += end + 3;
                        return Ok(Item::Set(set));
                    }
                    if posix_shaped(span.strip_prefix(b"^").unwrap_or(span)) {
                        return Err(Refusal::Perl(format!("POSIX class [:{written}:] unknown")));
                    }
                }
            }
            return Ok(Item::Char('['));
        }
        if byte != b'\\' {
            return Ok(Item::Char(char::from(byte)));
        }
        let Some(c) = self.peek() else {
            return Err(perl("Unmatched ["));
        };
        self.at += 1;
        Ok(match c {
            b'x' => Item::Char(self.hex()?),
            b'o' => Item::Char(self.braced_octal()?),
            b'0'..=b'7' => {
                self.at -= 1;
                Item::Char(self.octal())
            }
            b'N' => {
                if self.peek() != Some(b'{') {
                    return Err(perl(
                        "\\N in a character class must be a named character: \\N{...}",
                    ));
                }
                let mut chars = self.code_points()?;
                if chars.len() == 1 {
                    Item::Char(chars.remove(0))
                } else {
                    Item::Sequence(chars)
                }
            }
            b'p' | b'P' => Item::Set(self.property(c)?),
            b'e' => Item::Char('\x1B'),
            b'a' => Item::Char('\x07'),
            b'b' => Item::Char('\x08'),
            b't' => Item::Char('\t'),
            b'n' => Item::Char('\n'),
            b'r' => Item::Char('\r'),
            b'f' => Item::Char('\x0C'),
            b'c' => Item::Char(self.control()?),
            b'd' | b'D' | b'w' | b'W' | b's' | b'S' | b'h' | b'H' | b'v' | b'V' => {
                Item::Set(class_escape(c, flags, self.rule_set))
            }
            _ => Item::Char(char::from(c)),
        })
    }

    /// A group after its `(`: the group, a flag group, or `None` for a
    /// comment.
    fn group(&mut self, flags: &mut Flags) -> Result<Option<Tok>, Refusal> {
        if self.peek() == Some(b'*') {
            return Err(self.verb());
        }
        if self.peek() != Some(b'?') {
            let capturing = !flags.n;
            return self
                .group_body(
                    capturing,
                    if capturing { "(" } else { "(?:" }.to_string(),
                    *flags,
                )
                .map(Some);
        }
        self.at += 1;
        match self.peek() {
            Some(b'#') => self.comment().map(|()| None),
            Some(b'=' | b'!') => Err(Refusal::Unsupported(LOOK_AROUND.to_string())),
            Some(b'<') if matches!(self.peek_at(1), Some(b'=' | b'!')) => {
                Err(Refusal::Unsupported(LOOK_AROUND.to_string()))
            }
            Some(b'<') => self.named_group(b'>', *flags),
            Some(b'\'') => self.named_group(b'\'', *flags),
            Some(b'P') if self.peek_at(1) == Some(b'<') => {
                self.at += 1;
                self.named_group(b'>', *flags)
            }
            Some(b'P') => Err(Refusal::Unsupported(
                "a named back-reference or recursion (?P=...) is not supported".to_string(),
            )),
            Some(b'>') => Err(Refusal::Unsupported(
                "an atomic group (?>...) is not supported".to_string(),
            )),
            Some(b'|') => Err(Refusal::Unsupported(
                "a branch-reset group (?|...) is not supported".to_string(),
            )),
            Some(b'(') => Err(Refusal::Unsupported(
                "a conditional (?(...)...) is not supported".to_string(),
            )),
            Some(b'{') => Err(Refusal::Unsupported(
                "code in a pattern (?{...}) is not supported".to_string(),
            )),
            Some(b'R' | b'&' | b'+' | b'0'..=b'9') => Err(Refusal::Unsupported(
                "recursion into a group is not supported".to_string(),
            )),
            Some(b'-') if self.peek_at(1).is_some_and(|b| b.is_ascii_digit()) => Err(
                Refusal::Unsupported("recursion into a group is not supported".to_string()),
            ),
            _ => self.flag_group(flags),
        }
    }

    /// `(*...)` after the `(`: a backtracking control verb or an alpha
    /// assertion, which the crate lacks, or a name or shape Perl refuses,
    /// worded for a verb when the name has an uppercase letter or is empty
    /// and unknown, and for a `(*...)` construct otherwise. The name runs
    /// to the first `)` or `:` and an argument after the `:` to the first
    /// `)`; `(*MARK:` and `(*:` need one, and an alpha assertion needs the
    /// `:`.
    fn verb(&self) -> Refusal {
        const VERBS: [&str; 8] = [
            "ACCEPT", "FAIL", "F", "COMMIT", "PRUNE", "SKIP", "THEN", "MARK",
        ];
        const ASSERTIONS: [&str; 13] = [
            "pla",
            "plb",
            "nla",
            "nlb",
            "atomic",
            "sr",
            "asr",
            "positive_lookahead",
            "negative_lookahead",
            "positive_lookbehind",
            "negative_lookbehind",
            "script_run",
            "atomic_script_run",
        ];
        let rest = &self.pattern[self.at + 1..];
        let end = rest
            .iter()
            .position(|b| matches!(b, b')' | b':'))
            .unwrap_or(rest.len());
        let name = String::from_utf8_lossy(&rest[..end]);
        let upper = name.bytes().any(|b| b.is_ascii_uppercase());
        let colon = rest.get(end) == Some(&b':');
        let argument = if colon {
            rest[end + 1..].iter().position(|b| *b == b')')
        } else {
            None
        };
        if colon && argument.is_none() {
            return perl(if upper {
                "Unterminated verb pattern argument"
            } else {
                "Unterminated '(*...' argument"
            });
        }
        if !colon && end == rest.len() {
            return perl(if upper {
                "Unterminated verb pattern"
            } else {
                "Unterminated '(*...' construct"
            });
        }
        let mark = name == "MARK" || (name.is_empty() && colon);
        if mark && argument.unwrap_or(0) == 0 {
            return Refusal::Perl(format!("Verb pattern '{name}' has a mandatory argument"));
        }
        if mark || VERBS.contains(&name.as_ref()) {
            return Refusal::Unsupported(format!(
                "the backtracking control verb (*{name}...) is not supported"
            ));
        }
        if ASSERTIONS.contains(&name.as_ref()) {
            if !colon {
                return Refusal::Perl(format!("'(*{name}' requires a terminating ':'"));
            }
            return Refusal::Unsupported(format!(
                "the alpha assertion (*{name}:...) is not supported"
            ));
        }
        if upper || name.is_empty() {
            return Refusal::Perl(format!("Unknown verb pattern '{name}'"));
        }
        Refusal::Perl(format!("Unknown '(*...)' construct '{name}'"))
    }

    /// `(?<name>`, `(?'name'`, `(?P<name>`: a capturing group.
    fn named_group(&mut self, close: u8, flags: Flags) -> Result<Option<Tok>, Refusal> {
        self.at += 1;
        let Some(end) = self.pattern[self.at..].iter().position(|b| *b == close) else {
            return Err(perl("Sequence (?<... not terminated"));
        };
        let name = String::from_utf8_lossy(&self.pattern[self.at..self.at + end]).into_owned();
        self.at += end + 1;
        self.group_body(true, format!("(?P<{name}>"), flags)
            .map(Some)
    }

    /// `(?flags)` or `(?flags:...)` after the `?`.
    fn flag_group(&mut self, flags: &mut Flags) -> Result<Option<Tok>, Refusal> {
        let start = self.at;
        while self
            .peek()
            .is_some_and(|b| b.is_ascii_alphabetic() || b == b'^' || b == b'-')
        {
            self.at += 1;
        }
        let letters = &self.pattern[start..self.at];
        let stop = self.peek();
        if !matches!(stop, Some(b')' | b':')) {
            let mut shown = String::from_utf8_lossy(letters).into_owned();
            if let Some(stop) = stop {
                shown.push(char::from(stop));
            }
            return Err(Refusal::Perl(format!(
                "Sequence (?{shown}...) not recognized"
            )));
        }
        let mut inner = *flags;
        let text = apply_flags(letters, &mut inner)?;
        if inner.unicode || inner.ascii > 0 {
            self.charset_flag = true;
        }
        self.at += 1;
        if stop == Some(b')') {
            *flags = inner;
            return Ok(Some(Tok::Flags {
                text: if text.is_empty() {
                    String::new()
                } else {
                    format!("(?{text})")
                },
                flags: inner,
            }));
        }
        let open = if text.is_empty() {
            "(?:".to_string()
        } else {
            format!("(?{text}:")
        };
        self.group_body(false, open, inner).map(Some)
    }

    /// The branches of a group up to its `)`.
    fn group_body(&mut self, capturing: bool, open: String, flags: Flags) -> Result<Tok, Refusal> {
        if capturing {
            self.groups += 1;
        }
        self.depth += 1;
        if self.depth > NEST_LIMIT {
            return Err(Refusal::Unsupported(format!(
                "groups nested more than {NEST_LIMIT} deep are not supported"
            )));
        }
        let branches = self.branches(flags)?;
        if self.peek() != Some(b')') {
            return Err(perl("Unmatched ("));
        }
        self.at += 1;
        self.depth -= 1;
        Ok(Tok::Group(Group {
            capturing,
            open,
            branches,
            flags,
        }))
    }
}

/// Perl's flag letters applied to `flags`, with the crate's flag text
/// (`i-ms`, empty when the crate needs none). A leading `^` resets to the
/// defaults, and only positive letters follow it; `a`, `aa`, `d`, `u`, `n`,
/// `x`, `xx` and the no-op letters leave no crate text. The charset letters
/// `a`, `d`, `u` and `l` follow Perl's checks: `d` not after `^`, none after
/// the `-`, `a` at most twice and no other charset letter beside it, `d`, `u`
/// and `l` at most once and one of them; `l` passing them is the dialect's
/// refusal. A letter Perl does not know, a `^` past the first place, a `-`
/// after the `^` or a second `-` is Perl's refusal, which names the letters
/// up to it.
fn apply_flags(letters: &[u8], flags: &mut Flags) -> Result<String, Refusal> {
    let mut negating = false;
    let mut use_defaults = false;
    let unrecognized = |n: usize| {
        Refusal::Perl(format!(
            "Sequence (?{}...) not recognized",
            String::from_utf8_lossy(&letters[..=n])
        ))
    };
    let mut crate_flags: [(char, Option<bool>); 3] = [('i', None), ('m', None), ('s', None)];
    let set = |letter: char, on: bool, crate_flags: &mut [(char, Option<bool>); 3]| {
        if let Some(slot) = crate_flags.iter_mut().find(|(l, _)| *l == letter) {
            slot.1 = Some(on);
        }
    };
    // The charset letter read so far and, for `a`, how many times.
    let mut charset: Option<(u8, u8)> = None;
    let mut x_count = 0u8;
    for (n, &letter) in letters.iter().enumerate() {
        match letter {
            b'^' if n > 0 => return Err(unrecognized(n)),
            b'-' if negating || use_defaults => return Err(unrecognized(n)),
            b'^' => {
                use_defaults = true;
                flags.i = false;
                flags.x = 0;
                flags.ascii = 0;
                flags.unicode = false;
                flags.n = false;
                for slot in &mut crate_flags {
                    slot.1 = Some(false);
                }
            }
            b'-' => negating = true,
            b'i' => {
                flags.i = !negating;
                set('i', !negating, &mut crate_flags);
            }
            b'm' | b's' => set(char::from(letter), !negating, &mut crate_flags),
            b'x' => {
                x_count = x_count.saturating_add(1);
                flags.x = if negating { 0 } else { x_count.min(2) };
            }
            b'n' => flags.n = !negating,
            b'a' | b'd' | b'u' | b'l' => {
                let shown = char::from(letter);
                // Perl's order: a `d` after `^` is unrecognized; `a` and `d`
                // are tested against the `-` before an earlier charset letter,
                // `u` and `l` after it.
                if letter == b'd' && use_defaults {
                    return Err(unrecognized(n));
                }
                if negating && (matches!(letter, b'a' | b'd') || charset.is_none()) {
                    return Err(Refusal::Perl(format!(
                        "Regexp modifier \"{shown}\" may not appear after the \"-\""
                    )));
                }
                let count = match charset {
                    None => 1,
                    Some((b'a', 1)) if letter == b'a' => 2,
                    Some((b'a', _)) => {
                        return Err(perl("Regexp modifier \"a\" may appear a maximum of twice"))
                    }
                    Some((before, _)) if before == letter => {
                        return Err(Refusal::Perl(format!(
                            "Regexp modifier \"{shown}\" may not appear twice"
                        )))
                    }
                    Some((before, _)) => {
                        return Err(Refusal::Perl(format!(
                            "Regexp modifiers \"{}\" and \"{shown}\" are mutually exclusive",
                            char::from(before)
                        )))
                    }
                };
                charset = Some((letter, count));
                flags.ascii = if letter == b'a' { count } else { 0 };
                flags.unicode = letter == b'u';
            }
            b'p' | b'o' | b'g' | b'c' => {}
            _ => return Err(unrecognized(n)),
        }
    }
    if charset.is_some_and(|(letter, _)| letter == b'l') {
        return Err(Refusal::Unsupported(
            "the locale flag (?l) is not supported".to_string(),
        ));
    }
    let letters_set_to = |state: bool| -> String {
        crate_flags
            .iter()
            .filter(|(_, on)| *on == Some(state))
            .map(|(letter, _)| *letter)
            .collect()
    };
    let mut text = letters_set_to(true);
    let off = letters_set_to(false);
    if !off.is_empty() {
        text.push('-');
        text.push_str(&off);
    }
    Ok(text)
}

/// `\d \D \w \W \s \S \h \H \v \V` as the crate reads them under `flags`.
fn class_escape(c: u8, flags: Flags, rule_set: RuleSet) -> String {
    let ascii = rule_set.ascii_classes(flags);
    match c {
        b'd' if ascii => "[0-9]",
        b'D' if ascii => "[^0-9]",
        b'w' if ascii => "[0-9A-Za-z_]",
        b'W' if ascii => "[^0-9A-Za-z_]",
        b's' if ascii => r"[\t\n\x0B\x0C\r\x20]",
        b'S' if ascii => r"[^\t\n\x0B\x0C\r\x20]",
        b'd' => r"\d",
        b'D' => r"\D",
        b'w' => r"\w",
        b'W' => r"\W",
        b's' => r"\s",
        b'S' => r"\S",
        b'h' => r"[\t\x20\xA0]",
        b'H' => r"[^\t\x20\xA0]",
        b'v' => r"[\n\x0B\x0C\r\x85]",
        b'V' => r"[^\n\x0B\x0C\r\x85]",
        _ => unreachable!("a class escape"),
    }
    .to_string()
}

/// `\b` and `\B`: ASCII under `a`, `aa` and byte rules over the Latin-1
/// reading.
fn word_boundary(c: u8, flags: Flags, rule_set: RuleSet) -> String {
    if rule_set.ascii_classes(flags) {
        format!("(?-u:\\{})", char::from(c))
    } else {
        format!("\\{}", char::from(c))
    }
}

/// Whether Perl refuses `[:name:]` as an unknown POSIX class rather than
/// reading the `[` as a literal: a name of three to fourteen bytes with no
/// uppercase ASCII letter, no space or tab and no escaped `]`.
fn posix_shaped(name: &[u8]) -> bool {
    (3..=14).contains(&name.len())
        && !name
            .iter()
            .any(|&b| b.is_ascii_uppercase() || matches!(b, b' ' | b'\t'))
        && !name.windows(2).any(|w| w == b"\\]")
}

/// A POSIX class `[:name:]` as the crate reads it: the crate's ASCII class
/// over the bytes, under `a` or `aa` and under byte rules over the Latin-1
/// reading, the Latin-1-aware Unicode class otherwise; `None` for a name
/// Perl does not know. Case-insensitive, `lower` and `upper` are the cased
/// letters, which the ASCII `alpha` class and the `Cased` property hold.
fn posix_class(name: &str, negated: bool, flags: Flags, rule_set: RuleSet) -> Option<String> {
    const NAMES: [&str; 14] = [
        "alpha", "alnum", "ascii", "blank", "cntrl", "digit", "graph", "lower", "print", "punct",
        "space", "upper", "word", "xdigit",
    ];
    if !NAMES.contains(&name) {
        return None;
    }
    let cased = flags.i && matches!(name, "lower" | "upper");
    if rule_set.rules == Rules::Bytes || rule_set.ascii_classes(flags) {
        let caret = if negated { "^" } else { "" };
        let name = if cased { "alpha" } else { name };
        return Some(format!("[[:{caret}{name}:]]"));
    }
    let (positive, negative) = match name {
        "alpha" => (r"\p{Alphabetic}", r"\P{Alphabetic}"),
        "alnum" => (r"[\p{Alphabetic}\p{Nd}]", r"[^\p{Alphabetic}\p{Nd}]"),
        "ascii" => (r"[\x00-\x7F]", r"[^\x00-\x7F]"),
        "blank" => (r"[\t\p{Zs}]", r"[^\t\p{Zs}]"),
        "cntrl" => (r"\p{Cc}", r"\P{Cc}"),
        "digit" => (r"\p{Nd}", r"\P{Nd}"),
        "graph" => (r"[^\s\p{Cc}]", r"[\s\p{Cc}]"),
        "lower" | "upper" if cased => (r"\p{Cased}", r"\P{Cased}"),
        "lower" => (r"\p{Lowercase}", r"\P{Lowercase}"),
        "print" => (r"[^\p{Cc}]", r"\p{Cc}"),
        "punct" => (r"[\p{P}$+<=>^`|~]", r"[^\p{P}$+<=>^`|~]"),
        "space" => (r"\s", r"\S"),
        "upper" => (r"\p{Uppercase}", r"\P{Uppercase}"),
        "word" => (r"\w", r"\W"),
        "xdigit" => (r"[0-9A-Fa-f]", r"[^0-9A-Fa-f]"),
        _ => unreachable!("a POSIX class name"),
    };
    Some(if negated { negative } else { positive }.to_string())
}

/// A `\N{U+}` sequence of several characters where Perl reads the class
/// oddly, refused for the dialect.
fn sequence_in_class(place: &str) -> Refusal {
    Refusal::Unsupported(format!(
        "a \\N{{U+...}} escape of several characters {place} is not supported"
    ))
}

/// One item of a bracketed class as read.
#[derive(Debug)]
enum ClassItem {
    Char(char),
    /// A range from a character to a later one.
    Range(char, char),
    /// A class escape or a POSIX class, as the crate reads it.
    Set(String),
}

/// The bracketed class of `items` as the crate reads it. Under byte rules
/// over the Latin-1 reading a class with a character above ASCII is written
/// case-sensitive, with the other case of each ASCII letter and letter range
/// added when the match is case-insensitive, since the crate's own folding
/// would pair the Latin-1 letters Perl's byte rules leave uncased.
fn class_text(items: &[ClassItem], negated: bool, flags: Flags, rule_set: RuleSet) -> String {
    let rules = rule_set.rules;
    let above_ascii = items.iter().any(|item| match item {
        ClassItem::Char(ch) => !ch.is_ascii(),
        ClassItem::Range(_, hi) => !hi.is_ascii(),
        ClassItem::Set(_) => false,
    });
    let sensitive = rule_set.bytes_over_latin1(flags) && above_ascii;
    let fold_ascii = sensitive && flags.i;
    let other_case = |letter: char| char::from(u8::try_from(letter).expect("ASCII") ^ 0x20);
    let mut text = String::new();
    if sensitive {
        text.push_str("(?-i:");
    }
    text.push('[');
    if negated {
        text.push('^');
    }
    for item in items {
        match item {
            ClassItem::Char(ch) => {
                text.push_str(&class_char(*ch, rules));
                if fold_ascii && ch.is_ascii_alphabetic() {
                    text.push(other_case(*ch));
                }
            }
            ClassItem::Range(lo, hi) => {
                text.push_str(&class_char(*lo, rules));
                text.push('-');
                text.push_str(&class_char(*hi, rules));
                if fold_ascii {
                    let inside = |c: char| (*lo..=*hi).contains(&c);
                    for (first, last) in [('a', 'z'), ('A', 'Z')] {
                        let (from, to) = ((*lo).max(first), (*hi).min(last));
                        if from <= to && !(inside(other_case(from)) && inside(other_case(to))) {
                            text.push(other_case(from));
                            text.push('-');
                            text.push(other_case(to));
                        }
                    }
                }
            }
            ClassItem::Set(set) => text.push_str(set),
        }
    }
    text.push(']');
    if sensitive {
        text.push(')');
    }
    text
}

/// The characters the crate reads as syntax outside a class.
const META: &str = r"\.+*?()|[]{}^$#&-~";

/// A literal character as the crate reads it.
fn lit(ch: char, rules: Rules, out: &mut String) {
    if ch.is_ascii_alphanumeric() {
        out.push(ch);
    } else if META.contains(ch) {
        out.push('\\');
        out.push(ch);
    } else if ch.is_ascii_graphic() {
        out.push(ch);
    } else {
        push_hex(ch, rules, out);
    }
}

/// A character inside a bracketed class as the crate reads it.
fn class_char(ch: char, rules: Rules) -> String {
    let mut out = String::new();
    if ch.is_ascii_alphanumeric() {
        out.push(ch);
    } else {
        push_hex(ch, rules, &mut out);
    }
    out
}

/// The hex escape of `ch`: `\xHH` for a byte under byte rules, `\x{H..}`
/// under Unicode rules.
fn push_hex(ch: char, rules: Rules, out: &mut String) {
    let code = u32::from(ch);
    if rules == Rules::Bytes && code <= 0xFF {
        let _ = write!(out, "\\x{code:02X}");
    } else {
        let _ = write!(out, "\\x{{{code:X}}}");
    }
}

/// The full case fold of `ch` where it is longer than one character and
/// lies in Latin-1, `s` for `s`, `S` and the long s, and empty for every
/// other character.
fn full_fold(ch: char) -> &'static str {
    match ch {
        '\u{DF}' | '\u{1E9E}' => "ss",
        '\u{FB00}' => "ff",
        '\u{FB01}' => "fi",
        '\u{FB02}' => "fl",
        '\u{FB03}' => "ffi",
        '\u{FB04}' => "ffl",
        '\u{FB05}' | '\u{FB06}' => "st",
        's' | 'S' | '\u{17F}' => "s",
        _ => "",
    }
}

/// The fold set of `ch` as a key: its full fold where it has one, else its
/// lowercase, so that `[sS]`, `[\x{17F}s]` and `[\x{DF}\x{1E9E}]` are one
/// set each and `[xy]` is two.
fn fold_key(ch: char) -> String {
    let fold = full_fold(ch);
    if fold.is_empty() {
        ch.to_lowercase().collect()
    } else {
        fold.to_string()
    }
}

/// Whether `ch` lies above Latin-1 and folds to an ASCII letter, a fold
/// `aa` forbids: the long s and the Kelvin sign.
fn folds_to_ascii(ch: char) -> bool {
    matches!(ch, '\u{17F}' | '\u{212A}')
}

/// The longest run of `s` a fold regex is written for: the regex of a run
/// grows with the square of its length.
const LONGEST_RUN: usize = 64;

struct Emitter {
    rule_set: RuleSet,
    out: String,
    /// The longest run of `s` written so far.
    longest_run: usize,
}

/// A run of adjacent literals held back from the output, with the flags its
/// first character was read under: a literal under other fold flags ends it.
#[derive(Default)]
struct Run {
    chars: Vec<char>,
    flags: Flags,
}

impl Emitter {
    fn rules(&self) -> Rules {
        self.rule_set.rules
    }

    fn branches(&mut self, branches: &[Vec<Tok>], flags: Flags) {
        let mut flags = flags;
        for (n, branch) in branches.iter().enumerate() {
            if n > 0 {
                self.out.push('|');
            }
            let mut run = Run::default();
            self.items(branch, &mut flags, &mut run, None);
            self.flush(&mut run);
        }
    }

    /// The tokens of one branch onto `run` and the output; `joined` is
    /// `Some` for the branch of a group the run continues through
    /// ([`transparent`]), holding the `m` and `s` flags of the text the join
    /// drops, which every atom that reads them gets back.
    fn items(
        &mut self,
        items: &[Tok],
        flags: &mut Flags,
        run: &mut Run,
        mut joined: Option<Dropped>,
    ) {
        let mut i = 0;
        while i < items.len() {
            let tok = &items[i];
            let quant = match items.get(i + 1) {
                Some(Tok::Quant { text, never }) => Some((text.as_str(), *never)),
                _ => None,
            };
            i += if quant.is_some() { 2 } else { 1 };
            if quant.is_some_and(|(_, never)| never) {
                self.flush(run);
                self.out.push_str(match self.rules() {
                    Rules::Bytes => r"[^\x00-\xFF]",
                    Rules::Latin1 => r"[^\x00-\x{10FFFF}]",
                });
                continue;
            }
            let quant = quant.map(|(text, _)| text);
            match tok {
                Tok::Lit(ch) => match quant {
                    Some(quant) => {
                        self.flush(run);
                        self.quantified_lit(*ch, quant, *flags);
                    }
                    None => self.push(*ch, run, *flags),
                },
                Tok::Class(class) => {
                    if quant.is_none() && class.uniform && self.rule_set.folds(*flags) {
                        self.push(class.members[0], run, *flags);
                    } else {
                        self.flush(run);
                        self.class(class, quant, *flags);
                    }
                }
                Tok::Group(group) => {
                    if quant.is_none() && transparent(group, *flags) {
                        let mut inner = group.flags;
                        let mut dropped = joined.unwrap_or_default();
                        dropped.read(&group.open);
                        if let [branch] = group.branches.as_slice() {
                            self.items(branch, &mut inner, run, Some(dropped));
                        }
                    } else {
                        self.flush(run);
                        let dropped = joined.unwrap_or_default();
                        self.out.push_str(&dropped.open());
                        self.out.push_str(&group.open);
                        self.branches(&group.branches, group.flags);
                        self.out.push(')');
                        self.out.push_str(dropped.close());
                        self.out.push_str(quant.unwrap_or_default());
                    }
                }
                Tok::Flags { text, flags: after } => {
                    // A run continues through a flag group that keeps its
                    // fold state; the text then precedes the run, whose
                    // literals read neither `m` nor `s`.
                    if self.rule_set.fold_state(*after) != self.rule_set.fold_state(*flags) {
                        self.flush(run);
                    }
                    match &mut joined {
                        Some(dropped) => dropped.read(text),
                        None => self.out.push_str(text),
                    }
                    *flags = *after;
                }
                Tok::Other(text) => {
                    self.flush(run);
                    let dropped = joined.unwrap_or_default();
                    self.out.push_str(&dropped.open());
                    self.out.push_str(text);
                    self.out.push_str(dropped.close());
                    self.out.push_str(quant.unwrap_or_default());
                }
                Tok::Quant { .. } => unreachable!("a quantifier follows its atom"),
            }
        }
    }

    /// `ch` onto the run: its fold under full folding, itself otherwise. A
    /// run under other fold flags is written first.
    fn push(&mut self, ch: char, run: &mut Run, flags: Flags) {
        if !run.chars.is_empty()
            && self.rule_set.fold_state(run.flags) != self.rule_set.fold_state(flags)
        {
            self.flush(run);
        }
        if run.chars.is_empty() {
            run.flags = flags;
        }
        let fold = full_fold(ch);
        if self.rule_set.folds(flags) && !fold.is_empty() {
            run.chars.extend(fold.chars());
        } else {
            run.chars.push(ch);
        }
    }

    /// The run written out: as the regex of its fold under full folding,
    /// as literals otherwise.
    fn flush(&mut self, run: &mut Run) {
        if run.chars.is_empty() {
            return;
        }
        if self.rule_set.folds(run.flags) {
            self.fold_regex(&run.chars);
        } else {
            for &ch in &run.chars {
                self.literal(ch, run.flags);
            }
        }
        run.chars.clear();
    }

    /// A literal outside a fold run, kept case-sensitive where the crate
    /// would fold it and Perl does not: under `aa` a character above Latin-1
    /// that folds to ASCII, and under byte rules over the Latin-1 reading a
    /// character above ASCII.
    fn literal(&mut self, ch: char, flags: Flags) {
        let sensitive = flags.i
            && ((flags.ascii == 2 && folds_to_ascii(ch))
                || (!ch.is_ascii() && self.rule_set.bytes_over_latin1(flags)));
        if sensitive {
            self.out.push_str("(?-i:");
            lit(ch, self.rules(), &mut self.out);
            self.out.push(')');
        } else {
            lit(ch, self.rules(), &mut self.out);
        }
    }

    /// A quantified literal: under full folding a character with a fold
    /// string repeats the string or itself.
    fn quantified_lit(&mut self, ch: char, quant: &str, flags: Flags) {
        let fold = full_fold(ch);
        if self.rule_set.folds(flags) && fold.len() > 1 {
            self.out.push_str("(?:");
            self.fold_regex(&fold.chars().collect::<Vec<_>>());
            self.out.push(')');
        } else {
            self.literal(ch, flags);
        }
        self.out.push_str(quant);
    }

    /// A bracketed class: under full folding a listed member with a fold
    /// string also matches that string, and a `\N{U+}` sequence listed is an
    /// alternative, a run of literals of its own.
    fn class(&mut self, class: &Class, quant: Option<&str>, flags: Flags) {
        let member_folds: Vec<Vec<char>> = if self.rule_set.folds(flags) && !class.negated {
            class
                .members
                .iter()
                .map(|&member| full_fold(member))
                .filter(|fold| fold.len() > 1)
                .map(|fold| fold.chars().collect())
                .collect()
        } else {
            Vec::new()
        };
        if member_folds.is_empty() && class.sequences.is_empty() {
            self.out.push_str(&class.text);
        } else {
            self.out.push_str("(?:");
            self.out.push_str(&class.text);
            let alternatives = member_folds.iter().chain(class.sequences.iter());
            for (n, chars) in alternatives.enumerate() {
                if n > 0 || !class.text.is_empty() {
                    self.out.push('|');
                }
                let mut run = Run::default();
                for &ch in chars {
                    self.push(ch, &mut run, flags);
                }
                self.flush(&mut run);
            }
            self.out.push(')');
        }
        self.out.push_str(quant.unwrap_or_default());
    }

    /// The regex of a fold string: every character as itself, and every run
    /// of `k` letters `s` as the words over `s` and `ß` whose folds are `k`
    /// letters `s`, so that the fold of the match equals the fold of the run.
    /// A run past [`LONGEST_RUN`] is written plain and refused afterwards.
    fn fold_regex(&mut self, fold: &[char]) {
        let mut i = 0;
        while i < fold.len() {
            if fold[i] == 's' {
                let k = fold[i..].iter().take_while(|c| **c == 's').count();
                self.longest_run = self.longest_run.max(k);
                if k <= LONGEST_RUN {
                    s_run(k, &mut self.out);
                } else {
                    self.out.extend(std::iter::repeat_n('s', k));
                }
                i += k;
            } else {
                lit(fold[i], self.rules(), &mut self.out);
                i += 1;
            }
        }
    }
}

/// Whether a run of literals continues through `group`: it captures
/// nothing, has one branch or only empty ones, and neither its opening nor
/// a flag group inside changes case sensitivity, whose text the join would
/// drop. A charset change inside is admitted: it writes no text, and the
/// run splits where the fold flags change. A change of `m` or `s` inside is
/// admitted too: the join keeps it in a [`Dropped`] and writes it back
/// around every atom that reads it.
fn transparent(group: &Group, flags: Flags) -> bool {
    if group.capturing || group.flags.i != flags.i {
        return false;
    }
    match group.branches.as_slice() {
        [branch] => branch.iter().all(|tok| match tok {
            Tok::Flags { flags: after, .. } => after.i == group.flags.i,
            _ => true,
        }),
        branches => branches.iter().all(|b| b.is_empty()),
    }
}

/// The `m` and `s` flags of the group openings and flag groups a join
/// dropped, as the crate would read them: the last setting of each letter.
/// Only `.`, `^` and `$` read them, so they are written back around each
/// atom the join writes as text.
#[derive(Debug, Clone, Copy, Default)]
struct Dropped {
    m: Option<bool>,
    s: Option<bool>,
}

impl Dropped {
    /// The letters of a crate flag text, `(?s:`, `(?m-s)` or empty, read in.
    fn read(&mut self, text: &str) {
        let letters = text.trim_start_matches("(?").trim_end_matches([':', ')']);
        let (on, off) = letters.split_once('-').unwrap_or((letters, ""));
        for (letters, state) in [(on, true), (off, false)] {
            for letter in letters.chars() {
                match letter {
                    'm' => self.m = Some(state),
                    's' => self.s = Some(state),
                    _ => {}
                }
            }
        }
    }

    /// The opening of the group that gives the dropped flags back, `(?s:`;
    /// empty when none was dropped.
    fn open(self) -> String {
        let mut on = String::new();
        let mut off = String::new();
        for (flag, letter) in [(self.m, 'm'), (self.s, 's')] {
            match flag {
                Some(true) => on.push(letter),
                Some(false) => off.push(letter),
                None => {}
            }
        }
        match (on.is_empty(), off.is_empty()) {
            (true, true) => String::new(),
            (_, true) => format!("(?{on}:"),
            _ => format!("(?{on}-{off}:"),
        }
    }

    /// The closing of [`Dropped::open`]'s group; empty when it opened none.
    fn close(self) -> &'static str {
        if self.m.is_none() && self.s.is_none() {
            ""
        } else {
            ")"
        }
    }
}

/// The words over `s` (weight 1) and `ß` (weight 2) of weight `k`: a word
/// splits at the middle fold position, or a `ß` straddles it.
fn s_run(k: usize, out: &mut String) {
    match k {
        0 => {}
        1 => out.push('s'),
        2 => out.push_str(r"(?:ss|\x{DF})"),
        _ => {
            let a = k / 2;
            let b = k - a;
            out.push_str("(?:");
            s_run(a, out);
            s_run(b, out);
            out.push('|');
            s_run(a - 1, out);
            out.push_str(r"\x{DF}");
            s_run(b - 1, out);
            out.push(')');
        }
    }
}

/// `code` as a Unicode scalar value the crate accepts: a surrogate or a
/// value above U+10FFFF, which Perl admits in a pattern and no byte read as
/// Latin-1 ever matches, becomes the nearest scalar value below it, which no
/// such byte matches either and which keeps every range it bounds in order.
fn scalar(code: u32) -> char {
    let code = match code {
        0xD800..=0xDFFF => 0xD7FF,
        0x11_0000.. => 0x10_FFFF,
        _ => code,
    };
    char::from_u32(code).expect("a scalar value")
}

/// The code `value` names, or Perl's refusal of a value above its largest
/// code point, 0x7FFFFFFFFFFFFFFF. The refusal names the written `digits`
/// in their own `radix`; with `as_written` off, for a braced `\x{}` or
/// `\o{}` escape, a value below 0xFFFFFFFFFFFFFFFF, Perl's overflow mark,
/// is named in hex instead. A value above `u32` is kept as `u32::MAX`,
/// which [`scalar`] reads as a character no value holds.
fn code_point(value: u128, digits: &[u8], radix: u32, as_written: bool) -> Result<u32, Refusal> {
    const PERL_MAX: u128 = 0x7FFF_FFFF_FFFF_FFFF;
    if value > PERL_MAX {
        let (prefix, shown, max) = if !as_written && value < u128::from(u64::MAX) {
            ("0x", format!("{value:X}"), "7FFFFFFFFFFFFFFF")
        } else {
            let digits = String::from_utf8_lossy(digits).into_owned();
            match radix {
                16 => ("0x", digits, "7FFFFFFFFFFFFFFF"),
                _ => ("0", digits, "777777777777777777777"),
            }
        };
        return Err(Refusal::Perl(format!(
            "Use of code point {prefix}{shown} is not allowed; the permissible max is {prefix}{max}"
        )));
    }
    Ok(u32::try_from(value).unwrap_or(u32::MAX))
}

/// A braced `\x{}` or `\o{}` escape read as Perl reads it.
struct Braced {
    /// The value, saturated at `u128::MAX`.
    value: u128,
    /// The radix of the digits, 16 or 8.
    radix: u32,
    /// The span of the pattern the value was read from: the digits and the
    /// underscores among them, after the leading blanks.
    digits: std::ops::Range<usize>,
    /// The index after the `}`.
    after: usize,
    /// Whether the braces hold nothing, which `\o{}` refuses.
    empty: bool,
}

/// The braced escape whose `{` is `pattern[at]`: the content runs to the
/// first `}`; a space or tab at its start is skipped, an underscore before a
/// digit is skipped, the first other character that is no digit of `radix`
/// ends the value and the rest of the content is ignored, so no digit at all
/// is 0. `None` when `pattern[at]` is not `{` or no `}` follows, an escape
/// Perl refuses.
fn braced_number(pattern: &[u8], at: usize, radix: u32) -> Option<Braced> {
    if pattern.get(at) != Some(&b'{') {
        return None;
    }
    let close = at + pattern[at..].iter().position(|b| *b == b'}')?;
    let content = &pattern[at + 1..close];
    let start = content
        .iter()
        .take_while(|b| matches!(b, b' ' | b'\t'))
        .count();
    let mut i = start;
    let mut value: u128 = 0;
    while let Some(&byte) = content.get(i) {
        if let Some(digit) = char::from(byte).to_digit(radix) {
            value = value
                .saturating_mul(u128::from(radix))
                .saturating_add(u128::from(digit));
        } else if byte != b'_'
            || !content
                .get(i + 1)
                .is_some_and(|b| char::from(*b).is_digit(radix))
        {
            break;
        }
        i += 1;
    }
    Some(Braced {
        value,
        radix,
        digits: at + 1 + start..at + 1 + i,
        after: close + 1,
        empty: content.is_empty(),
    })
}

/// The value of the unbraced `\x` escape whose digits start at `pattern[at]`,
/// read as Perl reads it, with the index after them: the hex digits present,
/// two at most, so no digit at all is 0 and a third digit is a literal.
fn unbraced_hex(pattern: &[u8], at: usize) -> (u32, usize) {
    let digits = pattern[at..]
        .iter()
        .take(2)
        .take_while(|b| b.is_ascii_hexdigit())
        .count();
    let value = pattern[at..at + digits]
        .iter()
        .filter_map(|&b| char::from(b).to_digit(16))
        .fold(0, |value, digit| value * 16 + digit);
    (value, at + digits)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn crate_pattern(pattern: &[u8]) -> String {
        translate(pattern).unwrap().pattern
    }

    fn refused(pattern: &[u8]) -> Refusal {
        translate(pattern).unwrap_err()
    }

    /// The rules a pattern runs under: bytes unless a `\p{}` class, a wide `\x{}`,
    /// `\o{}` or octal escape, a `\N{}` escape or a `u`, `a` or `aa` flag appears,
    /// the flag outside a bracketed class and a comment; an `l` flag is refused.
    #[test]
    fn unicode_rules_follow_the_measured_triggers() {
        let latin1 = |p: &[u8]| translate(p).unwrap().rules == Rules::Latin1;
        assert!(!latin1(br"\\p{L}\x{00E9}"));
        assert!(latin1(br"a\x{0100}"));
        assert!(latin1(br"a\x{ 100}"));
        assert!(!latin1(br"a\x{+100}\x{00E9\p}"));
        assert!(latin1(br"r\o{400}?"));
        assert!(!latin1(br"r\o{377}"));
        assert!(latin1(br"r\400?"));
        assert!(!latin1(br"r\377"));
        assert!(latin1(br"\N{U+41}"));
        assert!(!latin1(br"\N"));
        assert!(latin1(b"(?u)ss"));
        assert!(latin1(b"(?a)ss"));
        assert!(latin1(b"(?aa:ss)"));
        assert!(!latin1(b"(?i)ss(?-i)"));
        assert_eq!(
            refused(b"(?l)ss"),
            Refusal::Unsupported("the locale flag (?l) is not supported".to_string())
        );
        assert!(!latin1(b"[(?u)]ss$"));
        assert!(!latin1(b"[[:alpha:](?u)]"));
        assert!(latin1(b"(?#[)(?u)ss"));
        assert!(!latin1(b"(?x)# (?u)\nss"));
        assert!(latin1(b"(?x)# c\n(?u)ss"));
        assert!(latin1(b"(?x)(?-x)# (?u)ss"));
        assert!(latin1(b"(?x:a)?# (?u)ss"));
        assert!(latin1(b"(?:(?x)a)# (?u)ss"));
        assert!(!latin1(b"(?x)(?:(?-x)a)?# (?u)"));
        assert!(latin1(b"(?:(?#c)(?x)a)?#(?u)ss"));
        assert!(!latin1(b"(?x)(?:(?#c)(?-x)a)?#(?u)"));
        assert_eq!(
            refused(b"(?#c"),
            Refusal::Perl("Sequence (?#... not terminated".to_string())
        );
        assert!(latin1(br"[\x{100}]"));
        assert!(!latin1(br"[\x{FF}](?i)"));
        assert!(!latin1(br"(?xx)[ ](?u)]"));
        assert_eq!(
            refused(b"[[:foo:]]"),
            Refusal::Perl("POSIX class [:foo:] unknown".to_string())
        );
        assert_eq!(
            refused(b"[[:^foo:]]"),
            Refusal::Perl("POSIX class [:^foo:] unknown".to_string())
        );
        assert_eq!(
            crate_pattern(b"[[:a:]][[:Alpha:]][[: alpha:]]"),
            r"[\x5B\x3Aa\x3A]\][\x5B\x3AAlpha\x3A]\][\x5B\x3A\x20alpha\x3A]\]"
        );
        assert_eq!(crate_pattern(b"[[:\xE9:]]"), r"[\x5B\x3A\xE9\x3A]\]");
        assert_eq!(
            crate_pattern(b"[[:abcdefghijklmno:]]"),
            r"[\x5B\x3Aabcdefghijklmno\x3A]\]"
        );
        assert_eq!(
            refused(b"[[:abcdefghijklmn:]]"),
            Refusal::Perl("POSIX class [:abcdefghijklmn:] unknown".to_string())
        );
        assert_eq!(
            refused(b"[[:^abcdefghijklmn:]]"),
            Refusal::Perl("POSIX class [:^abcdefghijklmn:] unknown".to_string())
        );
        assert_eq!(crate_pattern(b"[[:ab\t:]]"), r"[\x5B\x3Aab\x09\x3A]\]");
        assert_eq!(
            refused(b"[[:ab\n:]]"),
            Refusal::Perl("POSIX class [:ab\n:] unknown".to_string())
        );
        assert_eq!(crate_pattern(br"[[:a\]b:]]"), r"[\x5B\x3Aa\x5Db\x3A]\]");
        assert_eq!(
            refused(br"[[:a\d:]]"),
            Refusal::Perl(r"POSIX class [:a\d:] unknown".to_string())
        );
        // The groups are counted through a class that holds a parenthesis.
        assert_eq!(
            crate_pattern(br"[[:alpha:](x)]\1(y)"),
            r"[[[:alpha:]]\x28x\x29]\1(y)"
        );
        assert_eq!(
            refused(br"[[:alpha:](x)]\2(y)"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
    }

    /// A group nested past the nesting limit, a bound type, a backtracking
    /// verb, an escape past Perl's largest code point and a control escape of
    /// a byte above ASCII or of `{` are refused as Perl or the crate refuses
    /// them; Perl names the escape's own letter and the digits it read.
    #[test]
    fn deep_nesting_bound_types_verbs_and_wide_escapes_are_refused() {
        let nested = |n: usize| {
            let mut p = Vec::new();
            p.extend(std::iter::repeat_n(b"(?:".as_slice(), n).flatten());
            p.push(b's');
            p.extend(std::iter::repeat_n(b')', n));
            translate(&p)
        };
        assert!(nested(100).is_ok());
        assert_eq!(
            nested(101).unwrap_err(),
            Refusal::Unsupported("groups nested more than 100 deep are not supported".to_string())
        );
        assert!(nested(3000).is_err());
        assert_eq!(
            refused(br"\b{2}"),
            Refusal::Perl("'2' is an unknown bound type".to_string())
        );
        assert_eq!(
            refused(br"s\b{wb}"),
            Refusal::Unsupported(r"the bound type \b{wb} is not supported".to_string())
        );
        assert_eq!(
            refused(br"s\B{gcb}"),
            Refusal::Unsupported(r"the bound type \B{gcb} is not supported".to_string())
        );
        assert_eq!(
            refused(br"s\B{ gcb }"),
            Refusal::Perl("' gcb ' is an unknown bound type".to_string())
        );
        assert_eq!(
            refused(br"s\b{ }"),
            Refusal::Perl("' ' is an unknown bound type".to_string())
        );
        assert_eq!(
            refused(br"\b{wb"),
            Refusal::Perl(r"Missing right brace on \b{}".to_string())
        );
        assert_eq!(
            refused(br"\B{wb"),
            Refusal::Perl(r"Missing right brace on \B{}".to_string())
        );
        assert_eq!(refused(br"\b{}"), Refusal::Perl(r"Empty \b{}".to_string()));
        assert_eq!(refused(br"\B{}"), Refusal::Perl(r"Empty \B{}".to_string()));
        assert_eq!(
            refused(b"(*FAIL)s"),
            Refusal::Unsupported(
                "the backtracking control verb (*FAIL...) is not supported".to_string()
            )
        );
        assert_eq!(
            refused(b"(*:name)s"),
            Refusal::Unsupported(
                "the backtracking control verb (*...) is not supported".to_string()
            )
        );
        assert_eq!(
            refused(b"(*FAIL:)s"),
            Refusal::Unsupported(
                "the backtracking control verb (*FAIL...) is not supported".to_string()
            )
        );
        assert_eq!(
            refused(b"(*pla:s)"),
            Refusal::Unsupported("the alpha assertion (*pla:...) is not supported".to_string())
        );
        assert_eq!(
            refused(b"(*atomic_script_run:s)"),
            Refusal::Unsupported(
                "the alpha assertion (*atomic_script_run:...) is not supported".to_string()
            )
        );
        assert_eq!(
            refused(b"(*)s"),
            Refusal::Perl("Unknown verb pattern ''".to_string())
        );
        assert_eq!(
            refused(b"(*FOO)s"),
            Refusal::Perl("Unknown verb pattern 'FOO'".to_string())
        );
        assert_eq!(
            refused(b"(*xY)s"),
            Refusal::Perl("Unknown verb pattern 'xY'".to_string())
        );
        assert_eq!(
            refused(b"(*Foo:bar)s"),
            Refusal::Perl("Unknown verb pattern 'Foo'".to_string())
        );
        assert_eq!(
            refused(b"(*x)s"),
            Refusal::Perl("Unknown '(*...)' construct 'x'".to_string())
        );
        assert_eq!(
            refused(b"(*foo:bar)s"),
            Refusal::Perl("Unknown '(*...)' construct 'foo'".to_string())
        );
        assert_eq!(
            refused(b"(*MARK)s"),
            Refusal::Perl("Verb pattern 'MARK' has a mandatory argument".to_string())
        );
        assert_eq!(
            refused(b"(*MARK:)s"),
            Refusal::Perl("Verb pattern 'MARK' has a mandatory argument".to_string())
        );
        assert_eq!(
            refused(b"(*:)s"),
            Refusal::Perl("Verb pattern '' has a mandatory argument".to_string())
        );
        assert_eq!(
            refused(b"(*pla)"),
            Refusal::Perl("'(*pla' requires a terminating ':'".to_string())
        );
        assert_eq!(
            refused(b"(*asr)"),
            Refusal::Perl("'(*asr' requires a terminating ':'".to_string())
        );
        assert_eq!(
            refused(b"(*FAIL"),
            Refusal::Perl("Unterminated verb pattern".to_string())
        );
        assert_eq!(
            refused(b"(*FAIL:x"),
            Refusal::Perl("Unterminated verb pattern argument".to_string())
        );
        assert_eq!(
            refused(b"(*foo"),
            Refusal::Perl("Unterminated '(*...' construct".to_string())
        );
        assert_eq!(
            refused(b"(*"),
            Refusal::Perl("Unterminated '(*...' construct".to_string())
        );
        assert_eq!(
            refused(b"(*:"),
            Refusal::Perl("Unterminated '(*...' argument".to_string())
        );
        assert_eq!(
            refused(b"(*pla:s"),
            Refusal::Perl("Unterminated '(*...' argument".to_string())
        );
        // A reference to a missing group is Perl's refusal once the pattern
        // is read: before a construct the crate lacks, after a syntax error.
        // A `(` after the construct may open a group the first reading never
        // counted, so the reference is judged only when it exceeds them all.
        assert_eq!(
            refused(br"(r)\2(?=x)"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
        assert_eq!(
            refused(br"(r)\2(?=x)(s)"),
            Refusal::Unsupported(LOOK_AROUND.to_string())
        );
        assert_eq!(
            refused(br"(r)\2(?=x)(?:s)"),
            Refusal::Unsupported(LOOK_AROUND.to_string())
        );
        assert_eq!(
            refused(br"(r)\3(?=x)(s)"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
        assert_eq!(
            refused(br"(r)\2(?l)(s)"),
            Refusal::Unsupported("the locale flag (?l) is not supported".to_string())
        );
        assert_eq!(
            refused(br"(r)\2x)"),
            Refusal::Perl("Unmatched )".to_string())
        );
        assert_eq!(
            refused(br"(r)\2{01}"),
            Refusal::Perl("Invalid quantifier in {,}".to_string())
        );
        assert_eq!(
            refused(br"\c"),
            Refusal::Perl("Character following \"\\c\" must be printable ASCII".to_string())
        );
        assert_eq!(
            refused(br"\x{FFFFFFFFFFFFFFFFFFFFF}"),
            Refusal::Perl(
                "Use of code point 0xFFFFFFFFFFFFFFFFFFFFF is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF".to_string()
            )
        );
        assert_eq!(
            refused(br"\x{F_FFFFFFFFFFFFFFFFFFFFzz}"),
            Refusal::Perl(
                "Use of code point 0xF_FFFFFFFFFFFFFFFFFFFF is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF".to_string()
            )
        );
        assert_eq!(
            refused(br"\o{7777777777777777777777}"),
            Refusal::Perl(
                "Use of code point 07777777777777777777777 is not allowed; the permissible max is 0777777777777777777777".to_string()
            )
        );
        assert_eq!(
            refused(br"\o{1777777777777777777777}"),
            Refusal::Perl(
                "Use of code point 01777777777777777777777 is not allowed; the permissible max is 0777777777777777777777".to_string()
            )
        );
        // Below Perl's overflow mark the value itself is named, in hex for
        // either radix.
        assert_eq!(
            refused(br"\o{1000000000000000000001}"),
            Refusal::Perl(
                "Use of code point 0x8000000000000001 is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF".to_string()
            )
        );
        assert_eq!(
            refused(br"\x{00008000000000000000zz}"),
            Refusal::Perl(
                "Use of code point 0x8000000000000000 is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF".to_string()
            )
        );
        assert_eq!(crate_pattern(br"\x{7FFFFFFFFFFFFFFF}"), r"\x{10FFFF}");
        assert_eq!(crate_pattern(br"\o{777777777777777777777}"), r"\x{10FFFF}");
        assert_eq!(
            refused(b"\\c\xFF"),
            Refusal::Perl("Character following \"\\c\" must be printable ASCII".to_string())
        );
        assert_eq!(
            refused(br"\c{"),
            Refusal::Perl("Use \";\" instead of \"\\c{\"".to_string())
        );
        assert_eq!(crate_pattern(br"\c?\c@\c \c;"), r"\x7F\x00`\{");
    }

    /// Under byte rules a byte above ASCII is its `\xHH` escape, a backslash
    /// before one is dropped, and `\x{}` and `\x` read as Perl reads them.
    #[test]
    fn byte_patterns_keep_their_bytes() {
        assert_eq!(
            crate_pattern(b"r\\\xC3\xA9\\\\\xE9\\d"),
            r"r\xC3\xA9\\\xE9\d"
        );
        assert_eq!(crate_pattern(br"\xE$\xg\x"), r"\x0E$\x00g\x00");
        assert_eq!(crate_pattern(br"\x{00E9}\x{41}\x{zz}"), r"\xE9A\x00");
        assert_eq!(crate_pattern(br"\x{ 00E9 }\x{+00E9}"), r"\xE9\x00");
        assert_eq!(
            refused(br"\x{00E9"),
            Refusal::Perl("Missing right brace on \\x{}".to_string())
        );
        assert_eq!(crate_pattern(br"\e\a\cA\c?\t"), r"\x1B\x07\x01\x7F\x09");
        assert_eq!(crate_pattern(br"\y\Q\E\_\-"), r"yQE_\-");
        assert_eq!(crate_pattern(br"[\R\X\Z\y]"), r"[RXZy]");
        assert_eq!(crate_pattern(br"(?:|s)s(?:)s"), r"(?:|s)ss");
        let empty_branch = regex::bytes::Regex::new(&crate_pattern(br"^(?:s|)s$")).unwrap();
        assert!(empty_branch.is_match(b"s") && empty_branch.is_match(b"ss"));
        assert_eq!(
            crate_pattern(br"r\h\v\R"),
            r"r[\t\x20\xA0][\n\x0B\x0C\r\x85](?:\r\n|[\n\x0B\x0C\r\x85])"
        );
        assert_eq!(crate_pattern(br"\Aa\Z\K\G\N"), r"\Aa\z(?:)\A[^\n]");
    }

    /// Under `xx` a space or tab inside a bracketed class is skipped, before
    /// the `^`, before a leading `]` and after a range's `-` too; a newline
    /// is a member, and under `x` alone so is every blank.
    #[test]
    fn xx_skips_a_space_or_tab_inside_a_class() {
        assert_eq!(crate_pattern(b"(?xx)^[ ^a]$"), r"^[^a]$");
        assert_eq!(crate_pattern(b"(?xx)^[ ]]$"), r"^[\x5D]$");
        assert_eq!(crate_pattern(b"(?xx)^[^ ]]$"), r"^[^\x5D]$");
        assert_eq!(crate_pattern(b"(?xx)^[ ^\t]]$"), r"^[^\x5D]$");
        assert_eq!(crate_pattern(b"(?xx)^[a - z]$"), r"^[a-z]$");
        assert_eq!(crate_pattern(b"(?xx)^[a - ]$"), r"^[a\x2D]$");
        assert_eq!(crate_pattern(b"(?xx)^[a -]$"), r"^[a\x2D]$");
        assert_eq!(crate_pattern(b"(?xx)^[\ta]$"), r"^[a]$");
        assert_eq!(crate_pattern(b"(?xx)^[\na]$"), r"^[\x0Aa]$");
        assert_eq!(crate_pattern(b"(?xx)^[ \\ ]$"), r"^[\x20]$");
        assert_eq!(crate_pattern(b"(?x)^[ ^a]$"), r"^[\x20\x5Ea]$");
        assert_eq!(
            refused(b"(?xx)^[ ]$"),
            Refusal::Perl("Unmatched [".to_string())
        );
    }

    /// Under Unicode rules the bytes are Latin-1 characters, wide escapes are
    /// `\x{}` scalar values, and a surrogate or a value past U+10FFFF is the
    /// scalar value below it.
    #[test]
    fn latin1_patterns_read_bytes_as_characters() {
        assert_eq!(
            translate(b"\\\xE9\xE9\\x{ 0_100 }\\x{D800}\\x{FFFFFFFFF}").unwrap(),
            Translated {
                rules: Rules::Latin1,
                pattern: concat!(r"\x{E", "9}", r"\x{E", "9}", r"\x{100}\x{D7FF}\x{10FFFF}")
                    .to_string()
            }
        );
        assert_eq!(crate_pattern(br"\xE\p{L}\x"), r"\x{E}\p{L}\x{0}");
        assert_eq!(
            crate_pattern(br"r[\x{100}-\x{110000}]"),
            r"r[\x{100}-\x{10FFFF}]"
        );
    }

    /// Braces: a quantifier shape after an atom, a literal elsewhere; the
    /// refusals and the never-matching range.
    #[test]
    fn braces_quantify_or_stay_literal_as_perl_reads_them() {
        assert_eq!(
            crate_pattern(br"r{1,2}r{,3}r{ 1 , 2 }r{2}"),
            r"r{1,2}r{0,3}r{1,2}r{2}"
        );
        assert_eq!(
            crate_pattern(br"r{x}r{,}r{1,x}{"),
            r"r\{x\}r\{,\}r\{1,x\}\{"
        );
        assert_eq!(
            crate_pattern(br"(?:r){r+{r{1}{(?i){2}ss|{2}"),
            r"r\{r+\{r{1}(?i)\{\{2\}ss|\{2\}"
        );
        assert_eq!(crate_pattern(br"r{{1}"), r"r\{{1}");
        assert_eq!(crate_pattern(br"r{2,1}"), r"[^\x00-\xFF]");
        assert_eq!(crate_pattern(br"r{1,2}?r{1,2}+"), r"r{1,2}?r{1,2}");
        assert_eq!(
            refused(br"r{01}"),
            Refusal::Perl("Invalid quantifier in {,}".to_string())
        );
        assert_eq!(
            refused(br"r{0,02}"),
            Refusal::Perl("Invalid quantifier in {,}".to_string())
        );
        assert_eq!(
            refused(br"r{65535}"),
            Refusal::Perl("Quantifier in {,} bigger than 65534".to_string())
        );
        assert_eq!(
            refused(br"r{1}{2}"),
            Refusal::Perl("Nested quantifiers".to_string())
        );
        assert_eq!(
            refused(br"r**"),
            Refusal::Perl("Nested quantifiers".to_string())
        );
        assert_eq!(
            refused(br"r{1}??"),
            Refusal::Perl("Nested quantifiers".to_string())
        );
        assert_eq!(
            refused(br"a|*s"),
            Refusal::Perl("Quantifier follows nothing".to_string())
        );
        assert_eq!(crate_pattern(br"^{1}s(?:){2}"), r"^{1}s(?:){2}");
        assert_eq!(
            crate_pattern(br"a(?#c)+b(?#c){2}c+(?#c)+d(?#c)*?e(?#c)(?#d){"),
            r"a+b{2}c+d*?e\{"
        );
        assert_eq!(
            refused(br"a(?#c"),
            Refusal::Perl("Sequence (?#... not terminated".to_string())
        );
    }

    /// Octal: `\o{}`, `\0`, and the back-reference reading of `\1` to `\9`
    /// and of a longer number within the count of groups opened before it.
    #[test]
    fn octal_escapes_read_as_perl_does() {
        assert_eq!(
            crate_pattern(br"r\o{101}\o{ 101 }\o{1_01}\o{101x}\o{x}"),
            r"rAAAA\x00"
        );
        assert_eq!(
            crate_pattern(br"r\101\0\00\0001\0101\351\1010\18"),
            r"rA\x00\x00\x001\x081\xE9A0\x018"
        );
        assert_eq!(crate_pattern(br"(r)\10(r)\12"), r"(r)\x08(r)\x0A");
        assert_eq!(crate_pattern(br"(r)\1"), r"(r)\1");
        assert_eq!(crate_pattern(br"\1(r)"), r"\1(r)");
        let ten = "(r)".repeat(10);
        assert_eq!(
            crate_pattern(format!(r"\10{ten}\10\11").as_bytes()),
            format!(r"\x08{ten}\10\x09")
        );
        assert_eq!(
            refused(br"(r)\81"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
        assert_eq!(
            refused(br"r\1"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
        assert_eq!(
            refused(br"r\8"),
            Refusal::Perl("Reference to nonexistent group".to_string())
        );
        assert_eq!(refused(br"r\o{}"), Refusal::Perl("Empty \\o{}".to_string()));
        assert_eq!(
            refused(br"r\o{101"),
            Refusal::Perl("Missing right brace on \\o{}".to_string())
        );
        assert_eq!(
            refused(br"r\o101"),
            Refusal::Perl("Missing braces on \\o{}".to_string())
        );
        assert_eq!(
            crate_pattern(br"r[\101-\132][\0-\1][\o{101}]"),
            r"r[A-Z][\x00-\x01][A]"
        );
        assert_eq!(
            translate(br"\o{400}?\xe9").unwrap(),
            Translated {
                rules: Rules::Latin1,
                pattern: concat!(r"\x{100}?\x{E", "9}").to_string()
            }
        );
    }

    /// Full folding under Unicode rules: `ß` and `ss` run as the words whose
    /// fold is `ss`, the fold of a run is compared whole, and the joins and
    /// breaks of a run follow Perl.
    #[test]
    fn full_folds_rewrite_literal_runs() {
        let fold = |p: &[u8]| {
            let t = translate(p).unwrap();
            assert_eq!(t.rules, Rules::Latin1, "{}", String::from_utf8_lossy(p));
            t.pattern
        };
        let ss = r"(?:ss|\x{DF})";
        assert_eq!(fold(br"^ss\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(fold(br"^\x{DF}\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(fold(br"^\x{1E9E}\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(fold(b"^\xDF\\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(
            fold(br"^s\x{DF}s\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^sss\x{100}?$"),
            format!(r"^(?:s{ss}|\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(fold(br"^(?:s)s\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(fold(br"^s(?:)s\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(fold(br"^s(?:|)s\x{100}?$"), format!(r"^{ss}\x{{100}}?$"));
        assert_eq!(
            fold(br"^(?s:s)s(?m:s)[s]\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(fold(br"^(?s:s.)s\x{100}?$"), r"^s(?s:.)s\x{100}?$");
        assert_eq!(
            fold(br"^s(?s:(?:s))\x{100}?$"),
            format!(r"^{ss}\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?s:(?:s).)\x{100}?$"),
            format!(r"^{ss}(?s:.)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?m:(?:(?:s))$)(?-m:^)\x{100}?$"),
            format!(r"^{ss}(?m:$)(?-m:^)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?s:(?:s|x).)\x{100}?$"),
            r"^s(?s:(?:s|x))(?s:.)\x{100}?$"
        );
        assert_eq!(
            fold(br"^s(?:s(?s)s)s\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?:s(?m)(?i)s)s\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?:s(?-s)s)s\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?s:s(?x) s)s\x{100}?$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s(?:s(?s)s.)s\x{100}?$"),
            format!(r"^(?:s{ss}|\x{{DF}}s)(?s:.)s\x{{100}}?$")
        );
        let scoped = fold(br"^s(?:s(?:(?s)s).)s\x{100}?$");
        assert_eq!(scoped, format!(r"^(?:s{ss}|\x{{DF}}s).s\x{{100}}?$"));
        let scoped = regex::Regex::new(&format!("(?i){scoped}")).unwrap();
        assert!(scoped.is_match("\u{DF}s\ts") && !scoped.is_match("\u{DF}s\ns"));
        assert_eq!(
            fold(br"^s(?#c)s(?x) s\x{100}?$"),
            format!(r"^(?:s{ss}|\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^[s]s[sS\x{17F}]\x{100}?$"),
            format!(r"^(?:s{ss}|\x{{DF}}s)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^(ss)(s)(s)\x{100}?$"),
            format!(r"^({ss})(s)(s)\x{{100}}?$")
        );
        assert_eq!(fold(br"^(?:s|x)s\x{100}?$"), r"^(?:s|x)s\x{100}?$");
        assert_eq!(fold(br"^(?:s|)s\x{100}?$"), r"^(?:s|)s\x{100}?$");
        assert_eq!(fold(br"^(?:|s)s\x{100}?$"), r"^(?:|s)s\x{100}?$");
        assert_eq!(
            fold(br"^[xy][xX][s\x{17F}]s\x{100}?$"),
            format!(r"^[xy]x{ss}\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^(?:ss|zz)\x{100}?$"),
            format!(r"^(?:{ss}|zz)\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^s{2}ss{2}\x{DF}{2}\x{100}?$"),
            format!(r"^s{{2}}ss{{2}}(?:{ss}){{2}}\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^[\x{DF}]{2}[\x{DF}x][st]s\x{100}?$"),
            format!(r"^(?:[\x{{DF}}]|{ss}){{2}}(?:[\x{{DF}}x]|{ss})[st]s\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^[^\x{DF}][\x{DE}-\x{E0}]\x{100}?$"),
            r"^[^\x{DF}][\x{DE}-\x{E0}]\x{100}?$"
        );
        assert_eq!(
            fold(br"^(?-i)ss(?i)ss\x{100}?$"),
            format!(r"^(?-i)ss(?i){ss}\x{{100}}?$")
        );
        assert_eq!(
            fold(br"^(?-i:s)s(?aa)ss\x{17F}\x{100}?$"),
            r"^(?-i:s)sss(?-i:\x{17F})\x{100}?$"
        );
        assert_eq!(
            fold(br"^\x{FB01}f\x{FB01}\x{FB06}\x{100}?$"),
            r"^fiffist\x{100}?$"
        );
        assert_eq!(fold(br"(?^u)ss$"), "(?-ims)ss$");
        assert_eq!(fold(br"(?iu)ss$"), format!("(?i){ss}$"));
        assert_eq!(fold(br"(?a)^ss\w$"), format!("^{ss}[0-9A-Za-z_]$"));
        assert_eq!(fold(br"(?aa)^\xDF\x{DF}s$"), r"^\x{DF}\x{DF}s$");
        assert_eq!(
            fold(br"^(?u)ss(?n)(ss)$"),
            format!(r"^(?:{ss}{ss}|s\x{{DF}}s)$")
        );
        assert_eq!(fold(br"^s\b\Ks\x{100}?$"), r"^s\b(?:)s\x{100}?$");
        assert_eq!(
            fold(br"^[[:alpha:]][[:^lower:]](?a)[[:alpha:]]\x{100}?$"),
            r"^[\p{Alphabetic}][\P{Cased}][[[:alpha:]]]\x{100}?$"
        );
        assert_eq!(
            refused(b"(?l)ss"),
            Refusal::Unsupported("the locale flag (?l) is not supported".to_string())
        );
        assert_eq!(
            refused(b"(?e)ss"),
            Refusal::Perl("Sequence (?e...) not recognized".to_string())
        );
        assert_eq!(
            refused(b"s(?=s)"),
            Refusal::Unsupported(LOOK_AROUND.to_string())
        );
    }

    /// A run of `s` up to the bound compiles under full folding; one past it
    /// is refused naming the bound.
    #[test]
    fn a_run_past_the_bound_is_refused() {
        let run = |k: usize| {
            let mut pattern = br"\x{100}?".to_vec();
            pattern.extend(std::iter::repeat_n(b's', k));
            translate(&pattern)
        };
        let longest = run(LONGEST_RUN).unwrap();
        assert!(regex::Regex::new(&longest.pattern).is_ok());
        assert_eq!(
            run(LONGEST_RUN + 1).unwrap_err(),
            Refusal::Unsupported(
                "a run of more than 64 letters s under full case folding".to_string()
            )
        );
    }

    /// The words of weight `k` over `s` and `ß` are exactly the ones `s_run`
    /// matches: every split of `k` into ones and twos, nothing of another
    /// weight.
    #[test]
    fn s_runs_match_every_word_of_their_weight_alone() {
        let weight = |w: &str| {
            w.chars()
                .map(|c| if c == 's' { 1 } else { 2 })
                .sum::<usize>()
        };
        let mut words = vec![String::new()];
        let mut frontier = vec![String::new()];
        while !frontier.is_empty() {
            let mut next = Vec::new();
            for word in frontier {
                for c in ['s', '\u{DF}'] {
                    let mut longer = word.clone();
                    longer.push(c);
                    if weight(&longer) <= 12 {
                        words.push(longer.clone());
                        next.push(longer);
                    }
                }
            }
            frontier = next;
        }
        for k in 0..=12 {
            let mut out = String::from("^");
            s_run(k, &mut out);
            out.push('$');
            let re = regex::Regex::new(&out).unwrap();
            let matched: Vec<&String> = words.iter().filter(|w| re.is_match(w)).collect();
            assert!(matched.iter().all(|w| weight(w) == k), "k={k}: {matched:?}");
            assert_eq!(
                matched.len(),
                words.iter().filter(|w| weight(w) == k).count(),
                "k={k}"
            );
        }
    }

    /// A `u`, `a`, `aa` or `d` charset flag rules its own group or the rest
    /// of its group alone; outside it Perl's byte rules hold, written over the
    /// Latin-1 reading: ASCII classes, no case above ASCII and no fold to a
    /// sequence. A Unicode trigger anywhere puts the whole pattern under
    /// Unicode rules. A run of literals folds under the charset in force where
    /// it was read and splits where that changes, through adjacent groups.
    /// Case-insensitive `[[:lower:]]` and `[[:upper:]]` are the cased letters.
    /// The charset letters follow Perl's checks.
    #[test]
    fn charset_flags_rule_their_own_scope() {
        let latin1 = |p: &[u8]| translate(p).unwrap().rules == Rules::Latin1;
        assert!(!latin1(br"^\w\w$"));
        assert!(latin1(br"^(?u:\w)\w$"));
        assert!(latin1(br"^(?aa:ss)$"));
        assert_eq!(crate_pattern(br"^(?u:\w)\w$"), r"^\w[0-9A-Za-z_]$");
        assert_eq!(crate_pattern(br"^\w(?u:\w)$"), r"^[0-9A-Za-z_]\w$");
        assert_eq!(
            crate_pattern(br"^(?a:\w)\w$"),
            r"^[0-9A-Za-z_][0-9A-Za-z_]$"
        );
        assert_eq!(crate_pattern(br"^(?u:\w)\w\x{100}?$"), r"^\w\w\x{100}?$");
        assert_eq!(
            crate_pattern(br"^\w(?a:\w)\x{100}?$"),
            r"^\w[0-9A-Za-z_]\x{100}?$"
        );
        assert_eq!(crate_pattern(br"^\w(?u:\x{100}?)$"), r"^\w\x{100}?$");
        assert_eq!(crate_pattern(br"^(?d:\w)\w\x{100}?$"), r"^\w\w\x{100}?$");
        assert_eq!(crate_pattern(br"^\w(?u)\w$"), r"^[0-9A-Za-z_]\w$");
        assert_eq!(crate_pattern(br"^(?u)\w(?d)\w$"), r"^\w[0-9A-Za-z_]$");
        assert_eq!(
            crate_pattern(br"^(?u:(?a)\w)\w$"),
            r"^[0-9A-Za-z_][0-9A-Za-z_]$"
        );
        assert_eq!(
            crate_pattern(br"^(?^a:\w)\w$"),
            r"^(?-ims:[0-9A-Za-z_])[0-9A-Za-z_]$"
        );
        assert_eq!(
            crate_pattern(br"^(?u)(?^i)\w\xea$"),
            r"^(?i-ms)[0-9A-Za-z_](?-i:\x{EA})$"
        );
        assert_eq!(
            crate_pattern(br"^(?u:)\b\xea\B$"),
            r"^(?-u:\b)(?-i:\x{EA})(?-u:\B)$"
        );
        assert_eq!(
            crate_pattern(br"^(?u:)[[:alpha:]]\s$"),
            r"^[[[:alpha:]]][\t\n\x0B\x0C\r\x20]$"
        );
        assert_eq!(crate_pattern(br"^(?u:)[[:^word:]]$"), r"^[[[:^word:]]]$");
        // Case: a byte above ASCII has none, ASCII letters keep theirs.
        assert_eq!(crate_pattern(br"^(?u:a)\xea$"), r"^a(?-i:\x{EA})$");
        assert_eq!(crate_pattern(br"^(?u:)[^\xeaa]$"), r"^(?-i:[^\x{EA}aA])$");
        assert_eq!(crate_pattern(br"^(?u:)[s-\xea]$"), r"^(?-i:[s-\x{EA}S-Z])$");
        assert_eq!(crate_pattern(br"^(?u:)[A-\xc9]$"), r"^(?-i:[A-\x{C9}])$");
        assert_eq!(crate_pattern(br"^(?u:)[k-\xba]$"), r"^(?-i:[k-\x{BA}K-Z])$");
        assert_eq!(crate_pattern(br"^(?u:)[\xa0a]$"), r"^(?-i:[\x{A0}aA])$");
        assert_eq!(crate_pattern(br"^(?u:)[a-z]$"), r"^[a-z]$");
        assert_eq!(
            crate_pattern(br"^(?u:)(?-i)[\xea]$"),
            r"^(?-i)(?-i:[\x{EA}])$"
        );
        assert_eq!(
            crate_pattern(br"^(?u:)[[:lower:]][[:^upper:]]$"),
            r"^[[[:alpha:]]][[[:^alpha:]]]$"
        );
        assert_eq!(
            crate_pattern(br"^(?u:)(?-i)[[:lower:]]$"),
            r"^(?-i)[[[:lower:]]]$"
        );
        assert_eq!(
            crate_pattern(br"^[[:upper:]]\x{100}?$"),
            r"^[\p{Cased}]\x{100}?$"
        );
        assert_eq!(
            crate_pattern(br"^(?-i)[[:^upper:]]\x{100}?$"),
            r"^(?-i)[\P{Uppercase}]\x{100}?$"
        );
        // Folds: whole under one charset, split at a change, joined across
        // adjacent groups under the same one.
        let ss = r"(?:ss|\x{DF})";
        assert_eq!(
            crate_pattern(br"^(?u:ss)\xdf$"),
            format!("^{ss}(?-i:\\x{{DF}})$")
        );
        assert_eq!(
            crate_pattern(br"^\xdf(?u:ss)$"),
            format!("^(?-i:\\x{{DF}}){ss}$")
        );
        assert_eq!(crate_pattern(br"^(?aa:ss)\xdf$"), r"^ss(?-i:\x{DF})$");
        assert_eq!(
            crate_pattern(br"^(?a:ss)\xdf$"),
            format!("^{ss}(?-i:\\x{{DF}})$")
        );
        assert_eq!(crate_pattern(br"^(?u:s)s$"), "^ss$");
        assert_eq!(crate_pattern(br"^(?u:s)(?u:s)$"), format!("^{ss}$"));
        assert_eq!(crate_pattern(br"^(?u:s)(?aa:s)$"), "^ss$");
        assert_eq!(crate_pattern(br"^(?a:s)(?u:s)$"), format!("^{ss}$"));
        assert_eq!(
            crate_pattern(br"^(?u:s)s\x{100}?$"),
            format!("^{ss}\\x{{100}}?$")
        );
        assert_eq!(
            crate_pattern(br"^(?u:\xdf)(?u:\xdf)$"),
            format!("^(?:{ss}{ss}|s\\x{{DF}}s)$")
        );
        assert_eq!(crate_pattern(br"^(?u:)ss$"), "^ss$");
        assert_eq!(crate_pattern(br"^(?u:ss)(?n)(ss)$"), format!("^{ss}ss$"));
        assert_eq!(crate_pattern(br"^(?u:\xdf)$"), format!("^{ss}$"));
        assert_eq!(crate_pattern(br"^(?u:)\xdf+$"), r"^(?-i:\x{DF})+$");
        assert_eq!(crate_pattern(br"^(?u:x)|^ss$"), "^x|^ss$");
        // Perl's charset letter checks come before the locale refusal.
        let perl = |m: &str| Refusal::Perl(m.to_string());
        assert_eq!(
            refused(br"^(?-u:\w)$"),
            perl("Regexp modifier \"u\" may not appear after the \"-\"")
        );
        assert_eq!(
            refused(br"^(?-l:\w)$"),
            perl("Regexp modifier \"l\" may not appear after the \"-\"")
        );
        assert_eq!(
            refused(br"^(?a-d:\w)$"),
            perl("Regexp modifier \"d\" may not appear after the \"-\"")
        );
        assert_eq!(
            refused(br"^(?u-u:\w)$"),
            perl("Regexp modifier \"u\" may not appear twice")
        );
        assert_eq!(
            refused(br"^(?a-u:\w)$"),
            perl("Regexp modifier \"a\" may appear a maximum of twice")
        );
        assert_eq!(
            refused(br"^(?d-l:\w)$"),
            perl("Regexp modifiers \"d\" and \"l\" are mutually exclusive")
        );
        assert_eq!(
            refused(br"^(?^d:\w)$"),
            perl("Sequence (?^d...) not recognized")
        );
        assert_eq!(
            refused(br"^(?^id:\w)$"),
            perl("Sequence (?^id...) not recognized")
        );
        assert_eq!(
            refused(br"^(?^dd:\w)$"),
            perl("Sequence (?^d...) not recognized")
        );
        assert_eq!(
            refused(br"^(?i^:\w)$"),
            perl("Sequence (?i^...) not recognized")
        );
        assert_eq!(
            refused(br"^(?^-i:\w)$"),
            perl("Sequence (?^-...) not recognized")
        );
        assert_eq!(
            refused(br"^(?i-m-s:\w)$"),
            perl("Sequence (?i-m-...) not recognized")
        );
        assert_eq!(crate_pattern(br"^(?i-ms:\w)$"), r"^(?-ms:\w)$");
        let many_x = format!("^(?{}:a b)$", "x".repeat(300));
        assert_eq!(crate_pattern(many_x.as_bytes()), "^ab$");
        assert_eq!(
            refused(br"^(?ei:\w)$"),
            perl("Sequence (?e...) not recognized")
        );
        assert_eq!(
            refused(br"^(?d-u:\w)$"),
            perl("Regexp modifiers \"d\" and \"u\" are mutually exclusive")
        );
        assert_eq!(
            refused(br"^(?ua:\w)$"),
            perl("Regexp modifiers \"u\" and \"a\" are mutually exclusive")
        );
        assert_eq!(
            refused(br"^(?da:\w)$"),
            perl("Regexp modifiers \"d\" and \"a\" are mutually exclusive")
        );
        assert_eq!(
            refused(br"^(?lu:\w)$"),
            perl("Regexp modifiers \"l\" and \"u\" are mutually exclusive")
        );
        assert_eq!(
            refused(br"^(?uu:\w)$"),
            perl("Regexp modifier \"u\" may not appear twice")
        );
        assert_eq!(
            refused(br"^(?dd:\w)$"),
            perl("Regexp modifier \"d\" may not appear twice")
        );
        for letters in [&b"au"[..], b"ad", b"aaa", b"aua"] {
            let mut pattern = b"^(?".to_vec();
            pattern.extend_from_slice(letters);
            pattern.extend_from_slice(b":\\w)$");
            assert_eq!(
                refused(&pattern),
                perl("Regexp modifier \"a\" may appear a maximum of twice"),
                "{}",
                String::from_utf8_lossy(letters)
            );
        }
        assert_eq!(
            refused(br"^(?l:\w)$"),
            Refusal::Unsupported("the locale flag (?l) is not supported".to_string())
        );
    }

    /// `\N{U+...}`: hexadecimal numbers with single underscores between
    /// digits, separated by `.` into a sequence of characters that is one atom
    /// and one run of literals, an alternative inside a bracketed class; an
    /// empty or malformed number is Perl's refusal, a code point above Perl's
    /// largest is the escape refusal, and a sequence in a negated class, before
    /// a `-` or as a range end is the dialect's.
    #[test]
    fn named_escapes_read_code_point_sequences() {
        assert_eq!(crate_pattern(br"^r\N{U+41}$"), "^rA$");
        assert_eq!(crate_pattern(br"^r\N{U+4_1}$"), "^rA$");
        assert_eq!(crate_pattern(br"^r\N{U+0041.0042}$"), "^rAB$");
        assert_eq!(crate_pattern(br"^r\N{U+41.4_2}$"), "^rAB$");
        assert_eq!(crate_pattern(br"^r\N{U+41.42}{2}$"), "^r(?:AB){2}$");
        assert_eq!(crate_pattern(br"^r\N{U+41.42}(?#x)+?$"), "^r(?:AB)+?$");
        assert_eq!(crate_pattern(br"^r(?-i)\N{U+41.42}$"), "^r(?-i)AB$");
        assert_eq!(crate_pattern(br"^\N{U+73.73}$"), r"^(?:ss|\x{DF})$");
        assert_eq!(crate_pattern(br"^r[x\N{U+41.42}]$"), "^r(?:[x]|AB)$");
        assert_eq!(crate_pattern(br"^r[\N{U+41.42}]+$"), "^r(?:AB)+$");
        assert_eq!(crate_pattern(br"^r[\N{U+41.42}\N{U+42}]$"), "^r(?:[B]|AB)$");
        assert_eq!(crate_pattern(br"^r[\N{U+41}-\N{U+42}]$"), "^r[A-B]$");
        assert_eq!(crate_pattern(br"^r[-\N{U+41.42}]$"), r"^r(?:[\x{2D}]|AB)$");
        assert_eq!(
            refused(br"^r[\N{U+41.42}-"),
            Refusal::Perl("Unmatched [".to_string())
        );
        assert_eq!(
            refused(br"^r[^\N{U+41.42}-"),
            Refusal::Perl("Unmatched [".to_string())
        );
        assert_eq!(
            refused(br"^r[A-\N{U+41.42}"),
            Refusal::Perl("Unmatched [".to_string())
        );
        assert_eq!(crate_pattern(br"^r\N{U+41.42}\N{U+41}$"), "^rABA$");
        assert_eq!(crate_pattern(br"^r\N{U+0}\N{U+00000041}$"), r"^r\x{0}A$");
        assert_eq!(crate_pattern(br"^r\N{U+110000}$"), r"^r\x{10FFFF}$");
        let invalid = Refusal::Perl("Invalid hexadecimal number in \\N{U+...}".to_string());
        for pattern in [
            &br"^r\N{U+}$"[..],
            br"^r\N{U+zz}$",
            br"^r\N{U+41_}$",
            br"^r\N{U+_41}$",
            br"^r\N{U+4__1}$",
            br"^r\N{U+41.}$",
            br"^r\N{U+.41}$",
            br"^r\N{U+41..42}$",
            br"^r\N{U+-41}$",
            br"^r\N{U+x41}$",
            br"^r[\N{U+41.}]$",
        ] {
            assert_eq!(
                refused(pattern),
                invalid,
                "{}",
                String::from_utf8_lossy(pattern)
            );
        }
        assert_eq!(
            refused(br"^r\N{U+41"),
            Refusal::Perl("Missing right brace on \\N{}".to_string())
        );
        assert_eq!(
            refused(br"^r\N{U+800000000000000A}$"),
            Refusal::Perl(
                "Use of code point 0x800000000000000A is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF"
                    .to_string()
            )
        );
        assert_eq!(
            refused(br"^r\N{U+41.FFFFFFFFFFFFFFFFFFFF}$"),
            Refusal::Perl(
                "Use of code point 0xFFFFFFFFFFFFFFFFFFFF is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF"
                    .to_string()
            )
        );
        for pattern in [
            &br"^r[\N{U+41.42}-Z]$"[..],
            br"^r[\N{U+41.42}-]$",
            br"^r[\N{U+41.42}-Z-A]$",
        ] {
            assert_eq!(
                refused(pattern),
                Refusal::Unsupported(
                    "a \\N{U+...} escape of several characters before a `-` is not supported"
                        .to_string()
                ),
                "{}",
                String::from_utf8_lossy(pattern)
            );
        }
        for pattern in [
            &br"^r[\N{U+41.42}-B]$"[..],
            br"^r[\N{U+41.42}-B",
            br"^r[\N{U+41.42}-M]$",
        ] {
            assert_eq!(
                refused(pattern),
                Refusal::Perl("Invalid [] range".to_string()),
                "{}",
                String::from_utf8_lossy(pattern)
            );
        }
        assert_eq!(
            refused(br"^r[\N{U+41.42}-Z"),
            Refusal::Perl("Unmatched [".to_string())
        );
        assert_eq!(
            refused(br"^r\N{U+8_000_000_000_000_000}$"),
            Refusal::Perl(
                "Use of code point 0x8_000_000_000_000_000 is not allowed; the permissible max is 0x7FFFFFFFFFFFFFFF"
                    .to_string()
            )
        );
        assert_eq!(
            refused(br"^r[^\N{U+41.42}]$"),
            Refusal::Unsupported(
                "a \\N{U+...} escape of several characters in a negated class is not supported"
                    .to_string()
            )
        );
        assert_eq!(
            refused(br"^r[A-\N{U+41.42}]$"),
            Refusal::Unsupported(
                "a \\N{U+...} escape of several characters as the end of a range is not supported"
                    .to_string()
            )
        );
        assert_eq!(
            refused(br"^r[\N]$"),
            Refusal::Perl(
                "\\N in a character class must be a named character: \\N{...}".to_string()
            )
        );
    }
}
