// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The named items of a Rust source text and the lines each spans, which is how a record's
//! `implementation` table names code: `<path>::<item>`, the item a chain of names from the
//! file's top level through its modules, `impl` blocks and traits. The reader is lexical
//! and reads what `rust_items` of `scripts/ci/admission.py` reads;
//! `tests/fixtures/item_span_cases.jsonl` holds the two to one answer.

/// The words that may stand before an item's own keyword.
const MODIFIERS: [&str; 6] = ["pub", "async", "unsafe", "default", "auto", "safe"];
/// The words of an implemented type that are no part of its name.
const TYPE_QUALIFIERS: [&str; 5] = ["mut", "dyn", "const", "unsafe", "!"];
/// The prefixes of a string literal.
const STRING_PREFIXES: [&str; 5] = ["r", "br", "cr", "b", "c"];

/// What an item's header says it is.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    /// A module, a trait, an `impl` block or an `extern` block: its braces hold items.
    Container,
    /// A constant, a static, a type alias or an import: it ends at its semicolon, whatever
    /// braces its value holds.
    Statement,
    /// A function, a type, a macro: it ends at its closing brace or its semicolon.
    Block,
}

/// A word of an item's header with the count of generic brackets open around it.
type Word = (String, usize);

fn is_word_start(ch: char) -> bool {
    ch.is_ascii_alphabetic() || ch == '_'
}

fn is_word_char(ch: char) -> bool {
    ch.is_ascii_alphanumeric() || ch == '_'
}

fn starts_with(text: &[char], at: usize, prefix: &str) -> bool {
    prefix
        .chars()
        .enumerate()
        .all(|(offset, ch)| text.get(at + offset) == Some(&ch))
}

fn find(text: &[char], wanted: char, from: usize) -> Option<usize> {
    (from..text.len()).find(|&index| text[index] == wanted)
}

fn newlines(text: &[char], from: usize, to: usize) -> usize {
    text[from..to.min(text.len())]
        .iter()
        .filter(|&&ch| ch == '\n')
        .count()
}

/// The index after a string literal whose opening quote stands at `at`.
fn skip_string(text: &[char], at: usize) -> usize {
    let mut index = at + 1;
    while index < text.len() {
        match text[index] {
            '\\' => index += 2,
            '"' => return index + 1,
            _ => index += 1,
        }
    }
    text.len()
}

/// The index after a raw string literal whose hashes or opening quote stand at `at`, or
/// after the hashes alone where no quote follows them, as in a raw identifier.
fn skip_raw_string(text: &[char], at: usize) -> usize {
    let mut index = at;
    let mut hashes = 0;
    while text.get(index) == Some(&'#') {
        hashes += 1;
        index += 1;
    }
    if text.get(index) != Some(&'"') {
        return index;
    }
    let mut close = index + 1;
    while close < text.len() {
        if text[close] == '"' && (1..=hashes).all(|n| text.get(close + n) == Some(&'#')) {
            return close + 1 + hashes;
        }
        close += 1;
    }
    text.len()
}

/// The name of the type an `impl` header implements, from the words after `impl`: the
/// last word of the type's path, outside its generic brackets and before any `where`.
fn implemented_type(words: &[Word]) -> Option<String> {
    let mut outer: Vec<&str> = words
        .iter()
        .filter(|(_, angle)| *angle == 0)
        .map(|(word, _)| word.as_str())
        .collect();
    if let Some(place) = outer.iter().position(|word| *word == "where") {
        outer.truncate(place);
    }
    if let Some(place) = outer.iter().position(|word| *word == "for") {
        outer.drain(..=place);
    }
    outer
        .into_iter()
        .rfind(|word| !TYPE_QUALIFIERS.contains(word))
        .map(str::to_string)
}

/// What the words of an item's header make it, and its name where it has one.
fn classify(words: &[Word]) -> (Kind, Option<String>) {
    let word_at = |index: usize| words.get(index).map(|(word, _)| word.as_str());
    let named = |name: Option<&str>| name.filter(|name| *name != "!").map(str::to_string);
    let mut index = 0;
    while let Some(word) = word_at(index) {
        let following = word_at(index + 1);
        match word {
            _ if MODIFIERS.contains(&word) => {}
            "extern" => match following {
                Some("crate") => return (Kind::Statement, None),
                None => return (Kind::Container, None),
                Some(_) => {}
            },
            "const" if matches!(following, Some("fn" | "unsafe" | "async" | "extern")) => {}
            "const" => return (Kind::Statement, named(following)),
            "static" => {
                let name = if following == Some("mut") {
                    word_at(index + 2)
                } else {
                    following
                };
                return (Kind::Statement, named(name));
            }
            "type" => return (Kind::Statement, named(following)),
            "use" => return (Kind::Statement, None),
            "fn" | "struct" | "enum" | "union" => return (Kind::Block, named(following)),
            "mod" | "trait" => return (Kind::Container, named(following)),
            "impl" => return (Kind::Container, implemented_type(&words[index + 1..])),
            "macro_rules" => {
                let name = if following == Some("!") {
                    word_at(index + 2)
                } else {
                    None
                };
                return (Kind::Block, named(name));
            }
            _ => return (Kind::Block, None),
        }
        index += 1;
    }
    (Kind::Block, None)
}

/// An open module, trait or `impl` block.
struct Container {
    name: Option<String>,
    /// The brace depth of the items it holds.
    depth: usize,
    first: usize,
}

/// Open braces that hold no items.
struct Body {
    /// The brace depth outside them.
    depth: usize,
    /// The item they end, with its name and first line; `None` where the item goes on
    /// after them, to its semicolon.
    ends: Option<(Option<String>, usize)>,
}

/// Every named item of a Rust source text that its top level reaches through modules,
/// traits and `impl` blocks, as its name, its first line and its last, in the order the
/// items end. A name is the chain from the top level joined by `::`: `function`,
/// `Type::method` for a method of any `impl` block of the type, `Trait::method`,
/// `tests::function`. An item's first line is that of its first doc comment or attribute,
/// its last that of its closing brace or semicolon. Lines are counted by line feeds.
pub fn items(text: &str) -> Vec<(String, usize, usize)> {
    let text: Vec<char> = text.chars().collect();
    let n = text.len();
    let mut found = Vec::new();
    let mut containers = vec![Container {
        name: None,
        depth: 0,
        first: 1,
    }];
    let mut body: Option<Body> = None;
    let mut depth = 0usize;
    let mut line = 1usize;
    // The header of the item being read: its first line, its words outside parentheses,
    // and the parentheses and generic brackets open in it.
    let mut start: Option<usize> = None;
    let mut words: Vec<Word> = Vec::new();
    let mut parens = 0usize;
    let mut angle = 0usize;
    let mut i = 0usize;

    fn record(
        found: &mut Vec<(String, usize, usize)>,
        containers: &[Container],
        name: Option<String>,
        first: usize,
        last: usize,
    ) {
        let Some(name) = name.filter(|name| name != "_") else {
            return;
        };
        let mut chain = Vec::with_capacity(containers.len());
        for container in &containers[1..] {
            match &container.name {
                Some(outer) => chain.push(outer.as_str()),
                None => return,
            }
        }
        chain.push(name.as_str());
        found.push((chain.join("::"), first, last));
    }

    while i < n {
        let ch = text[i];
        if ch == '\n' {
            line += 1;
            i += 1;
        } else if starts_with(&text, i, "//") {
            let outer_doc = starts_with(&text, i, "///") && !starts_with(&text, i, "////");
            if outer_doc && body.is_none() && start.is_none() {
                start = Some(line);
            }
            i = find(&text, '\n', i).unwrap_or(n);
        } else if starts_with(&text, i, "/*") {
            let outer_doc = starts_with(&text, i, "/**")
                && !starts_with(&text, i, "/***")
                && !starts_with(&text, i, "/**/");
            if outer_doc && body.is_none() && start.is_none() {
                start = Some(line);
            }
            let mut level = 1;
            let mut j = i + 2;
            while j < n && level > 0 {
                if starts_with(&text, j, "/*") {
                    level += 1;
                    j += 2;
                } else if starts_with(&text, j, "*/") {
                    level -= 1;
                    j += 2;
                } else {
                    j += 1;
                }
            }
            line += newlines(&text, i, j);
            i = j.min(n);
        } else if ch == '"' {
            let j = skip_string(&text, i);
            line += newlines(&text, i, j);
            i = j;
        } else if ch == '\'' {
            if text.get(i + 1) == Some(&'\\') {
                // An escaped character: the quote that closes it stands after the escape.
                let j = find(&text, '\'', i + 3).map_or(n, |close| close + 1);
                line += newlines(&text, i, j);
                i = j;
            } else if text.get(i + 2) == Some(&'\'') {
                line += newlines(&text, i, i + 3);
                i += 3;
            } else {
                // A lifetime or a label: its name is no word of the header.
                i += 1;
                while i < n && is_word_char(text[i]) {
                    i += 1;
                }
            }
        } else if is_word_start(ch) {
            let mut j = i;
            while j < n && is_word_char(text[j]) {
                j += 1;
            }
            let word: String = text[i..j].iter().collect();
            if STRING_PREFIXES.contains(&word.as_str()) && matches!(text.get(j), Some('"' | '#')) {
                let after = if word == "b" || word == "c" {
                    if text[j] == '"' {
                        skip_string(&text, j)
                    } else {
                        j
                    }
                } else {
                    skip_raw_string(&text, j)
                };
                if after > j {
                    line += newlines(&text, i, after);
                    i = after;
                    continue;
                }
            }
            if body.is_none() {
                start.get_or_insert(line);
                if parens == 0 {
                    words.push((word, angle));
                }
            }
            i = j;
        } else if ch == '#' && (starts_with(&text, i, "#[") || starts_with(&text, i, "#![")) {
            let inner = text[i + 1] == '!';
            let mut j = i + if inner { 3 } else { 2 };
            let mut level = 1;
            while j < n && level > 0 {
                match text[j] {
                    '"' => {
                        j = skip_string(&text, j);
                        continue;
                    }
                    '[' => level += 1,
                    ']' => level -= 1,
                    _ => {}
                }
                j += 1;
            }
            if !inner && body.is_none() && start.is_none() {
                start = Some(line);
            }
            line += newlines(&text, i, j);
            i = j.min(n);
        } else if ch == '{' {
            if body.is_some() {
                // Braces inside braces that hold no items.
            } else if parens > 0 || angle > 0 {
                // A block inside the header's parentheses or generic brackets is a constant
                // expression, not the item's body.
                body = Some(Body { depth, ends: None });
            } else {
                let (kind, name) = classify(&words);
                let first = start.unwrap_or(line);
                match kind {
                    Kind::Container => {
                        containers.push(Container {
                            name,
                            depth: depth + 1,
                            first,
                        });
                        start = None;
                        words.clear();
                        parens = 0;
                        angle = 0;
                    }
                    Kind::Statement => body = Some(Body { depth, ends: None }),
                    Kind::Block => {
                        body = Some(Body {
                            depth,
                            ends: Some((name, first)),
                        });
                    }
                }
            }
            depth += 1;
            i += 1;
        } else if ch == '}' {
            i += 1;
            if depth == 0 {
                continue;
            }
            depth -= 1;
            let mut ended = false;
            if let Some(open) = &body {
                if depth == open.depth {
                    if let Some((name, first)) = body.take().and_then(|open| open.ends) {
                        record(&mut found, &containers, name, first, line);
                        ended = true;
                    }
                }
            } else if containers.len() > 1 && depth < containers[containers.len() - 1].depth {
                let closed = containers.pop().expect("an open container");
                record(&mut found, &containers, closed.name, closed.first, line);
                ended = true;
            }
            if ended {
                start = None;
                words.clear();
                parens = 0;
                angle = 0;
            }
        } else if ch == ';' {
            if body.is_none() && parens == 0 {
                if let Some(first) = start {
                    let (_, name) = classify(&words);
                    record(&mut found, &containers, name, first, line);
                }
                start = None;
                words.clear();
                angle = 0;
            }
            i += 1;
        } else {
            if body.is_none() {
                match ch {
                    '(' | '[' => parens += 1,
                    ')' | ']' => parens = parens.saturating_sub(1),
                    '<' if parens == 0 => angle += 1,
                    '>' if parens == 0
                        && angle > 0
                        && i > 0
                        && !matches!(text[i - 1], '-' | '=') =>
                    {
                        angle -= 1;
                    }
                    '!' if parens == 0 => words.push(("!".to_string(), angle)),
                    _ => {}
                }
            }
            i += 1;
        }
    }
    found
}

/// The first and the last line of the first item of the name in a Rust source text, or
/// `None` where the text defines no such item.
pub fn item_span(text: &str, item: &str) -> Option<(usize, usize)> {
    items(text)
        .into_iter()
        .find(|(name, _, _)| name == item)
        .map(|(_, first, last)| (first, last))
}

#[cfg(test)]
#[path = "item_tests.rs"]
mod tests;
