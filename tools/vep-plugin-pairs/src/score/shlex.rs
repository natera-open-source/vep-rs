// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The words of a command line as CPython's `shlex.split` reads them: split at spaces,
//! tabs, carriage returns and line feeds; single quotes keep everything between them;
//! double quotes keep everything but a backslash before a double quote or a backslash;
//! outside quotes a backslash keeps the character after it; no comments.

/// What the reader is inside of.
#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    Between,
    Word,
    Quote(char),
    /// After a backslash, with the state to return to.
    Escape(Return),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Return {
    Word,
    Quote(char),
}

fn is_separator(character: char) -> bool {
    matches!(character, ' ' | '\t' | '\r' | '\n')
}

/// The words of `line`, or the reason `shlex.split` raises: `No closing quotation` or
/// `No escaped character`.
pub fn split(line: &str) -> Result<Vec<String>, &'static str> {
    let mut words = Vec::new();
    let mut word = String::new();
    let mut state = State::Between;
    for character in line.chars() {
        state = match state {
            State::Between if is_separator(character) => State::Between,
            State::Between | State::Word => match character {
                separator if is_separator(separator) => {
                    words.push(std::mem::take(&mut word));
                    State::Between
                }
                '\\' => State::Escape(Return::Word),
                '\'' | '"' => State::Quote(character),
                other => {
                    word.push(other);
                    State::Word
                }
            },
            State::Quote(quote) => {
                if character == quote {
                    // The word stands even when the quotes held nothing.
                    State::Word
                } else if character == '\\' && quote == '"' {
                    State::Escape(Return::Quote(quote))
                } else {
                    word.push(character);
                    State::Quote(quote)
                }
            }
            State::Escape(back) => {
                // Inside double quotes a backslash stays unless it stands before the quote
                // or another backslash.
                if let Return::Quote(quote) = back {
                    if character != '\\' && character != quote {
                        word.push('\\');
                    }
                }
                word.push(character);
                match back {
                    Return::Word => State::Word,
                    Return::Quote(quote) => State::Quote(quote),
                }
            }
        };
    }
    match state {
        State::Quote(_) => Err("No closing quotation"),
        State::Escape(_) => Err("No escaped character"),
        State::Word => {
            words.push(word);
            Ok(words)
        }
        State::Between => Ok(words),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(line: &str) -> Vec<String> {
        split(line).unwrap()
    }

    #[test]
    fn command_line_is_split_as_a_shell_reads_it() {
        assert_eq!(
            words("vep --plugin 'Planted,a file' --offline"),
            ["vep", "--plugin", "Planted,a file", "--offline"]
        );
        assert_eq!(words("  a\tb\r\nc  "), ["a", "b", "c"]);
        assert_eq!(words(""), Vec::<String>::new());
        assert_eq!(words("   "), Vec::<String>::new());
        assert_eq!(words("a '' b"), ["a", "", "b"]);
        assert_eq!(words("''"), [""]);
        assert_eq!(words("a\"b c\"d"), ["ab cd"]);
        assert_eq!(words("a#b # c"), ["a#b", "#", "c"]);
        assert_eq!(words("a\\ b"), ["a b"]);
        assert_eq!(words("\"a\\\"b\\\\c\\d\""), ["a\"b\\c\\d"]);
        assert_eq!(words("'a\\'"), ["a\\"]);
        assert_eq!(words("a\\\nb"), ["a\nb"]);
        assert_eq!(words("é\u{a0}x"), ["é\u{a0}x"]);
    }

    #[test]
    fn unbalanced_quote_and_trailing_backslash_are_refused() {
        assert_eq!(split("vep 'a file"), Err("No closing quotation"));
        assert_eq!(split("vep \"a"), Err("No closing quotation"));
        assert_eq!(split("vep \\"), Err("No escaped character"));
        assert_eq!(split("\"a\\"), Err("No escaped character"));
    }
}
