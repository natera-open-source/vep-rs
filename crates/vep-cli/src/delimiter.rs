// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `--delimiter` value as Ensembl VEP's column-based parsers read it: a
//! regular expression Perl's `split` splits every data line on (ensembl-io
//! ColumnBasedParser.pm:250). [`Delimiter`] splits the patterns format
//! detection chooses and plain text itself; any other value is compiled here
//! and split as Perl's `split` splits on a pattern.

use std::ops::Range;
use std::sync::Arc;

use vep_io::input::{Delimiter, PatternSplit};

/// The delimiter a `--delimiter` value names ([`Delimiter::from_flag`]), a
/// value with a regular-expression character compiled as the pattern Perl's
/// `split` reads it as. A value the pattern syntax refuses, or one that could
/// match part of a UTF-8 character, is an error naming the value and the fault.
/// So is a pattern that can match the empty string: Perl's `split` then takes,
/// at each position, the first alternative whose match is not empty before it
/// moves on, a backtracking order a search for the leftmost match cannot
/// reproduce, so `x*`, `,?`, the bare `|` and `\b` are refused rather than
/// split differently; a pattern that always consumes a character splits as
/// Perl splits, since both searches take the leftmost match with the same
/// alternative and quantifier preferences.
pub(crate) fn from_flag(value: &str) -> Result<Delimiter, String> {
    Delimiter::from_flag(value, |pattern| {
        let split = RegexSplit::new(pattern).map_err(|fault| {
            format!("--delimiter {pattern:?} is not a regular expression this run can split on: {fault}")
        })?;
        if let Some(assertion) = assertion_in(pattern) {
            return Err(format!(
                "--delimiter {pattern:?} holds the zero-width assertion `{assertion}`, which can match the empty string; a delimiter pattern has to match at least one character"
            ));
        }
        if split.text.is_match("") {
            return Err(format!(
                "--delimiter {pattern:?} matches the empty string; a delimiter pattern has to match at least one character"
            ));
        }
        Ok(Arc::new(split) as Arc<dyn PatternSplit>)
    })
}

/// The first zero-width assertion `pattern` holds outside a character class:
/// `^`, `$`, `\A`, `\z`, `\b` (its `\b{...}` forms included), `\B`, `\<` or
/// `\>`; `None` without one. An assertion is the one construct that matches
/// the empty string at a position without matching the empty text; a pattern
/// without one matches the empty string anywhere exactly when it matches the
/// empty text, which [`from_flag`] tests. An escape stands for its character,
/// a class runs from `[` to its `]`, nested classes and a leading literal `]`
/// included, and a comment under the `x` flag is read like the pattern.
fn assertion_in(pattern: &str) -> Option<String> {
    let mut chars = pattern.chars().peekable();
    let mut class_depth = 0usize;
    while let Some(c) = chars.next() {
        match c {
            '\\' => match chars.next() {
                Some(escaped @ ('A' | 'z' | 'b' | 'B' | '<' | '>')) if class_depth == 0 => {
                    return Some(format!("\\{escaped}"));
                }
                _ => {}
            },
            '[' => {
                class_depth += 1;
                chars.next_if_eq(&'^');
                chars.next_if_eq(&']');
            }
            ']' if class_depth > 0 => class_depth -= 1,
            '^' | '$' if class_depth == 0 => return Some(c.to_string()),
            _ => {}
        }
    }
    None
}

/// A `--delimiter` pattern compiled over text and over bytes, Unicode off in
/// both, as Perl reads a line it opened without a decoding layer: `.`, a class
/// and `\s`, `\w`, `\d` match single bytes, a literal outside ASCII matches
/// its UTF-8 bytes, and a class or property over characters outside ASCII is
/// refused. The text form refuses a pattern that could match part of a UTF-8
/// character (`.`, a negated class, a byte escape above ASCII), so its matches
/// end on character boundaries; the byte form splits a line that is not UTF-8
/// as Perl splits the bytes it read.
struct RegexSplit {
    text: regex::Regex,
    bytes: regex::bytes::Regex,
}

impl RegexSplit {
    fn new(pattern: &str) -> Result<Self, String> {
        let text = regex::RegexBuilder::new(pattern)
            .unicode(false)
            .build()
            .map_err(fault)?;
        let bytes = regex::bytes::RegexBuilder::new(pattern)
            .unicode(false)
            .build()
            .map_err(fault)?;
        Ok(Self { text, bytes })
    }
}

/// The one line of a `regex` error that names the fault: the last line of a
/// syntax report, which spans several.
fn fault(error: regex::Error) -> String {
    match &error {
        regex::Error::Syntax(report) => report
            .lines()
            .last()
            .unwrap_or_default()
            .trim_start_matches("error: ")
            .to_string(),
        other => other.to_string(),
    }
}

impl PatternSplit for RegexSplit {
    fn is_match(&self, line: &[u8]) -> bool {
        self.bytes.is_match(line)
    }

    fn split_str<'a>(&self, line: &'a str) -> Vec<&'a str> {
        let spans = perl_split(
            line.len(),
            |at| {
                let found = self.text.captures_at(line, at)?;
                groups(found.len(), |i| found.get(i).map(|m| m.range()))
            },
            |at| line[at..].chars().next().map_or(1, char::len_utf8),
        );
        spans.into_iter().map(|span| &line[span]).collect()
    }

    fn split_bytes<'a>(&self, line: &'a [u8]) -> Vec<&'a [u8]> {
        let spans = perl_split(
            line.len(),
            |at| {
                let found = self.bytes.captures_at(line, at)?;
                groups(found.len(), |i| found.get(i).map(|m| m.range()))
            },
            |_| 1,
        );
        spans.into_iter().map(|span| &line[span]).collect()
    }
}

/// A match as [`perl_split`] reads it: the whole match's span, then each
/// capture group's, `None` for a group that took no part.
type Groups = (Range<usize>, Vec<Option<Range<usize>>>);

/// The spans of a match with `count` groups, group 0 the whole match.
fn groups(count: usize, span: impl Fn(usize) -> Option<Range<usize>>) -> Option<Groups> {
    Some((span(0)?, (1..count).map(span).collect()))
}

/// Perl's `split /pattern/, $line` (`pp_split`) over a line of `len` units:
/// `captures_at(at)` is the leftmost match starting at or after `at`, and
/// `step(at)` the width of the character at `at`. The fields are the text
/// before each match and the text after the last, a leading empty field kept
/// when a match opens the line and trailing empty fields dropped. A match must
/// end past the field's start (Perl's `minend` of 1), so an empty match there
/// is passed over for the next match one character on; a pattern
/// [`from_flag`] accepts yields none, and the rule keeps the loop advancing
/// for any other. Each capture group of a match is a field of its own, an
/// unmatched group an empty one (Perl's `undef`).
fn perl_split(
    len: usize,
    captures_at: impl Fn(usize) -> Option<Groups>,
    step: impl Fn(usize) -> usize,
) -> Vec<Range<usize>> {
    let mut fields: Vec<Range<usize>> = Vec::new();
    let mut start = 0;
    while start < len {
        let mut found = captures_at(start);
        if found.as_ref().is_some_and(|(whole, _)| whole.end == start) {
            found = captures_at(start + step(start));
        }
        let Some((whole, groups)) = found else {
            break;
        };
        fields.push(start..whole.start);
        fields.extend(groups.into_iter().map(|group| group.unwrap_or(0..0)));
        start = whole.end;
    }
    if start < len {
        fields.push(start..len);
    }
    while fields.last().is_some_and(|field| field.is_empty()) {
        fields.pop();
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pattern(text: &str) -> Delimiter {
        let delimiter = from_flag(text).unwrap();
        assert!(matches!(delimiter, Delimiter::Pattern(_)), "{text:?}");
        delimiter
    }

    /// The fields of `line` under `pattern`, which its text and its bytes
    /// split alike.
    fn split(pattern_text: &str, line: &str) -> Vec<String> {
        let delimiter = pattern(pattern_text);
        let fields = delimiter.split(line);
        let as_bytes: Vec<&[u8]> = fields.iter().map(|field| field.as_bytes()).collect();
        assert_eq!(
            delimiter.split_bytes(line.as_bytes()),
            as_bytes,
            "{pattern_text:?} {line:?}"
        );
        fields.into_iter().map(str::to_string).collect()
    }

    /// Perl's `split` on each pattern, every expected list its own: the fields
    /// between matches; a leading empty field for a match that opens the line
    /// and none for a trailing one; a capture group is a field, an unmatched
    /// one empty; `\s` is ASCII whitespace, so a no-break space is text.
    #[test]
    fn pattern_splits_as_perl_split() {
        assert_eq!(
            split("[ ,]", "21,100 100,A/G 1"),
            ["21", "100", "100", "A/G", "1"]
        );
        assert_eq!(split("[ ,]", ",a,,b, ,"), ["", "a", "", "b"]);
        assert_eq!(split("[ ,]", " 21 100"), ["", "21", "100"]);
        assert_eq!(split("\\s*,\\s*", "a , b,c"), ["a", "b", "c"]);
        assert_eq!(split("\\s+,", "a\u{a0} ,b"), ["a\u{a0}", "b"]);
        assert_eq!(split("\\|", "a|b|"), ["a", "b"]);
        assert!(split("[ ,]", "").is_empty());
        assert_eq!(split("(,)", "a,b"), ["a", ",", "b"]);
        assert_eq!(split("(,)|(;)", "a,b;c"), ["a", ",", "", "b", "", ";", "c"]);
        assert_eq!(split("[ ,]", "caf\u{e9},b"), ["café", "b"]);
    }

    /// A line of bytes that is not UTF-8 splits on the pattern's byte matches,
    /// each field the bytes between them.
    #[test]
    fn pattern_splits_bytes_that_are_not_utf8() {
        let delimiter = pattern("[ ,]");
        assert_eq!(
            delimiter.split_bytes(b"21,caf\xe9 C"),
            [&b"21"[..], b"caf\xe9", b"C"]
        );
    }

    /// Format detection finds the pattern in a line it matches and, in a line
    /// it does not, chooses as it does for any other delimiter.
    #[test]
    fn pattern_is_found_by_format_detection() {
        let delimiter = pattern("[,;]");
        assert_eq!(delimiter.chosen_for("a;b"), delimiter);
        assert_eq!(delimiter.chosen_for("ab"), delimiter);
        assert_eq!(delimiter.chosen_for("a  b"), Delimiter::SpaceRun);
    }

    /// A value the pattern syntax refuses is an error naming the value and the
    /// fault; so is a pattern that could match part of a UTF-8 character (`.`,
    /// a negated class, a byte escape above ASCII) or a class or property over
    /// characters outside ASCII, which Perl would read as bytes (a literal
    /// outside ASCII matches its bytes on both sides and is accepted); so is one
    /// that can match the empty string: one the empty text matches (`x*`, `,?`,
    /// the bare `|`, `(|a)`, `\s*`) or one holding an assertion (`\b`, `^`, `$`,
    /// `\A`, `\z`, `\B`), hidden behind a consuming alternative or not, the
    /// error naming the assertion. A `^` or `$` inside a character class, an
    /// escaped one and an escaped backslash before `b` are characters, so those
    /// patterns are accepted. A value with no regular-expression character is
    /// literal text.
    #[test]
    fn invalid_pattern_is_an_error_naming_it() {
        let error = from_flag("[").unwrap_err();
        assert!(
            error.starts_with(
                "--delimiter \"[\" is not a regular expression this run can split on: "
            ),
            "{error}"
        );
        assert!(error.ends_with("unclosed character class"), "{error}");
        assert!(from_flag("(?=x)").unwrap_err().contains("look-around"));
        for invalid_utf8 in ["(?-u:\\xC3)", ".", "[^,]", "\\xe9", "\\W"] {
            let error = from_flag(invalid_utf8).unwrap_err();
            assert!(error.contains("invalid UTF-8"), "{invalid_utf8:?}: {error}");
        }
        for unicode in ["[\u{e9}]", "\\p{Zs}"] {
            let error = from_flag(unicode).unwrap_err();
            assert!(error.contains("Unicode"), "{unicode:?}: {error}");
        }
        for empty in ["x*", ",?", "|", "(|a)", "\\s*"] {
            assert_eq!(
                from_flag(empty).unwrap_err(),
                format!("--delimiter {empty:?} matches the empty string; a delimiter pattern has to match at least one character")
            );
        }
        for (assertion, named) in [
            ("\\b", "\\b"),
            ("(?:\\b|,)", "\\b"),
            ("[,;]|^\\b|c", "^"),
            ("^", "^"),
            (",$", "$"),
            ("(?m)$", "$"),
            ("\\B", "\\B"),
            ("\\A,", "\\A"),
            ("[,]\\z", "\\z"),
        ] {
            assert_eq!(
                from_flag(assertion).unwrap_err(),
                format!("--delimiter {assertion:?} holds the zero-width assertion `{named}`, which can match the empty string; a delimiter pattern has to match at least one character")
            );
        }
        for consuming in [
            "[]^$]",
            "[[:space:]^]+",
            "\\^|\\$",
            "a\\\\b",
            "\u{e9}|\u{df}",
        ] {
            pattern(consuming);
        }
        assert_eq!(from_flag("::"), Ok(Delimiter::Literal("::".to_string())));
        assert_eq!(from_flag("\\t"), Ok(Delimiter::Tab));
    }
}
