// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! File names read from standard input, in the four ways a shell pipeline hands a list of
//! names to a digest program.

/// How the names stand on standard input.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Form {
    /// One a line, whole, as `while IFS= read -r f` reads them.
    Lines,
    /// One a line, the blanks at either end dropped, as `while read -r f` reads them.
    Read,
    /// Each ended by a zero byte, as `xargs -0` reads them.
    Nul,
    /// Split at blanks and newlines, quotes and backslashes read, as `xargs` reads them.
    Words,
}

impl Form {
    pub fn parse(word: &str) -> Option<Form> {
        match word {
            "lines" => Some(Form::Lines),
            "read" => Some(Form::Read),
            "nul" => Some(Form::Nul),
            "words" => Some(Form::Words),
            _ => None,
        }
    }

    /// Whether the names reach the digest program through `xargs`, which runs it once on
    /// no name at all and leaves with status 123 where a name fails.
    pub fn through_xargs(self) -> bool {
        matches!(self, Form::Nul | Form::Words)
    }
}

fn is_blank(byte: u8) -> bool {
    byte == b' ' || byte == b'\t'
}

/// The lines of `text` that end in a newline; a shell `read` loop does not run its body
/// on a last line without one.
fn ended_lines(text: &[u8]) -> impl Iterator<Item = &[u8]> {
    let whole = memchr::memrchr(b'\n', text).map_or(0, |last| last + 1);
    text[..whole]
        .split_inclusive(|byte| *byte == b'\n')
        .map(|line| &line[..line.len() - 1])
}

/// The names `xargs` makes of `text`: blanks and newlines separate them, a quote of
/// either kind holds what stands up to its match on the same line, and a backslash holds
/// the byte after it.
fn words(text: &[u8]) -> Result<Vec<Vec<u8>>, String> {
    let mut names = Vec::new();
    let mut name: Option<Vec<u8>> = None;
    let mut bytes = text.iter().copied();
    while let Some(byte) = bytes.next() {
        match byte {
            b'\n' | b' ' | b'\t' => names.extend(name.take()),
            b'\'' | b'"' => {
                let held = name.get_or_insert_with(Vec::new);
                loop {
                    match bytes.next() {
                        Some(next) if next == byte => break,
                        Some(b'\n') | None => {
                            let kind = if byte == b'"' { "double" } else { "single" };
                            return Err(format!(
                                "unmatched {kind} quote in the names on standard input"
                            ));
                        }
                        Some(next) => held.push(next),
                    }
                }
            }
            b'\\' => {
                let held = name.get_or_insert_with(Vec::new);
                if let Some(next) = bytes.next() {
                    held.push(next);
                }
            }
            other => name.get_or_insert_with(Vec::new).push(other),
        }
    }
    names.extend(name);
    Ok(names)
}

/// The names `text` holds under `form`, in its order.
pub fn parse(text: &[u8], form: Form) -> Result<Vec<Vec<u8>>, String> {
    Ok(match form {
        Form::Lines => ended_lines(text).map(<[u8]>::to_vec).collect(),
        Form::Read => ended_lines(text)
            .map(|line| {
                let start = line.iter().take_while(|byte| is_blank(**byte)).count();
                let end = line.len()
                    - line[start..]
                        .iter()
                        .rev()
                        .take_while(|byte| is_blank(**byte))
                        .count();
                line[start..end].to_vec()
            })
            .collect(),
        Form::Nul => {
            let mut names: Vec<Vec<u8>> =
                text.split(|byte| *byte == 0).map(<[u8]>::to_vec).collect();
            // What stands after the last zero byte is a name only where it holds a byte.
            if names.last().is_some_and(Vec::is_empty) {
                names.pop();
            }
            names
        }
        Form::Words => words(text)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(text: &[u8], form: Form) -> Vec<String> {
        parse(text, form)
            .unwrap()
            .into_iter()
            .map(|name| String::from_utf8(name).unwrap())
            .collect()
    }

    // while IFS= read -r f; do printf '[%s]\n' "$f"; done
    #[test]
    fn whole_lines_keep_their_blanks_and_a_last_line_without_a_newline_is_not_read() {
        assert_eq!(names(b"x  \n\ty\n\nlast", Form::Lines), ["x  ", "\ty", ""]);
        assert_eq!(names(b"", Form::Lines), [""; 0]);
        assert_eq!(names(b"only", Form::Lines), [""; 0]);
    }

    // while read -r f; do printf '[%s]\n' "$f"; done
    #[test]
    fn read_lines_lose_the_blanks_at_either_end() {
        assert_eq!(
            names(b"x  \n\ty\n  z z  \n \t \nlast", Form::Read),
            ["x", "y", "z z", ""]
        );
    }

    // xargs -0 sh -c 'for a; do printf "[%s]\n" "$a"; done' _
    #[test]
    fn zero_ended_names_keep_an_empty_name_and_a_last_name_without_its_zero() {
        assert_eq!(names(b"a\0\0b\0", Form::Nul), ["a", "", "b"]);
        assert_eq!(names(b"a b\0c\nd", Form::Nul), ["a b", "c\nd"]);
        assert_eq!(names(b"", Form::Nul), [""; 0]);
    }

    // xargs sh -c 'for a; do printf "[%s]\n" "$a"; done' _
    #[test]
    fn words_are_split_at_blanks_and_read_quotes_and_backslashes() {
        assert_eq!(
            names(b"'q u o'\nb\\ c\n\"d  e\" f\n\n  lead\nlast", Form::Words),
            ["q u o", "b c", "d  e", "f", "lead", "last"]
        );
        assert_eq!(names(b"a''b \"\"\n", Form::Words), ["ab", ""]);
        assert!(parse(b"it's\n", Form::Words).is_err());
        assert!(parse(b"\"open\nclosed\"\n", Form::Words).is_err());
    }
}
