// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Regular expressions as CPython's `re` module reads and runs them on a `str`: the syntax,
//! the messages of a pattern it refuses, the match `re.search` finds and the text `re.sub`,
//! `re.split` and `re.finditer` give. The patterns of the reproduction records are written
//! for that module, so their meaning here is its meaning.

mod compile;
mod parse;
mod vm;

pub use parse::PatternError;

use compile::Program;
use vm::Vm;

#[derive(Debug)]
pub struct Regex {
    program: Program,
}

/// One match: the spans of the whole match and of each group, in characters, with the text
/// they are spans of.
#[derive(Clone, Debug)]
pub struct Match {
    spans: Vec<Option<(usize, usize)>>,
    names: Vec<(String, usize)>,
    text: std::rc::Rc<[char]>,
}

impl Match {
    pub fn start(&self) -> usize {
        self.spans[0].expect("the whole match").0
    }

    pub fn end(&self) -> usize {
        self.spans[0].expect("the whole match").1
    }

    /// `m.group(index)`: None for a group that took no part in the match.
    pub fn group(&self, index: usize) -> Option<String> {
        let (start, end) = (*self.spans.get(index)?)?;
        Some(self.text[start..end].iter().collect())
    }

    /// `m.group(name)`.
    pub fn name(&self, name: &str) -> Option<String> {
        let index = self.names.iter().find(|(n, _)| n == name)?.1;
        self.group(index)
    }
}

impl Regex {
    /// `re.compile(pattern)`.
    pub fn new(pattern: &str) -> Result<Regex, PatternError> {
        let parsed = parse::parse(pattern)?;
        let program = compile::compile(&parsed)?;
        Ok(Regex { program })
    }

    fn found(&self, vm: &Vm, text: &std::rc::Rc<[char]>) -> Match {
        Match {
            spans: (0..self.program.groups).map(|g| vm.group(g)).collect(),
            names: self.program.names.clone(),
            text: text.clone(),
        }
    }

    /// The first match at or after `from`; with `must_advance`, one that is not the empty
    /// match at `from` itself.
    fn search_from(&self, vm: &mut Vm, from: usize, must_advance: bool, len: usize) -> bool {
        (from..=len).any(|start| {
            if let Some(first) = self.program.first {
                if !vm.starts_in(first, start) {
                    return false;
                }
            }
            vm.match_at(start, must_advance && start == from, false)
                .is_some()
        })
    }

    /// `re.search(pattern, text) is not None`.
    pub fn is_match(&self, text: &str) -> bool {
        let chars: Vec<char> = text.chars().collect();
        let mut vm = Vm::new(&self.program, &chars);
        self.search_from(&mut vm, 0, false, chars.len())
    }

    /// `re.search(pattern, text)`.
    pub fn search(&self, text: &str) -> Option<Match> {
        let chars: std::rc::Rc<[char]> = text.chars().collect();
        let mut vm = Vm::new(&self.program, &chars);
        self.search_from(&mut vm, 0, false, chars.len())
            .then(|| self.found(&vm, &chars))
    }

    /// `re.match(pattern, text)`: a match that starts at the text's start.
    pub fn match_start(&self, text: &str) -> Option<Match> {
        let chars: std::rc::Rc<[char]> = text.chars().collect();
        let mut vm = Vm::new(&self.program, &chars);
        vm.match_at(0, false, false)
            .is_some()
            .then(|| self.found(&vm, &chars))
    }

    /// `re.fullmatch(pattern, text) is not None`.
    pub fn is_full_match(&self, text: &str) -> bool {
        let chars: Vec<char> = text.chars().collect();
        let mut vm = Vm::new(&self.program, &chars);
        vm.match_at(0, false, true).is_some()
    }

    /// `list(re.finditer(pattern, text))`: the matches left to right, none overlapping, an
    /// empty match allowed right after a match that is not empty.
    pub fn find_iter(&self, text: &str) -> Vec<Match> {
        let chars: std::rc::Rc<[char]> = text.chars().collect();
        let mut vm = Vm::new(&self.program, &chars);
        let mut matches = Vec::new();
        let mut from = 0;
        let mut must_advance = false;
        while self.search_from(&mut vm, from, must_advance, chars.len()) {
            let found = self.found(&vm, &chars);
            must_advance = found.start() == found.end();
            from = found.end();
            matches.push(found);
        }
        matches
    }

    /// `re.findall(pattern, text)` for a pattern without groups: the matched texts.
    pub fn find_all(&self, text: &str) -> Vec<String> {
        self.find_iter(text)
            .iter()
            .map(|m| m.group(0).expect("the whole match"))
            .collect()
    }

    /// `re.sub(pattern, replacement, text)` for a replacement that holds no backslash, so
    /// it stands for itself.
    pub fn sub(&self, replacement: &str, text: &str) -> String {
        let matches = self.find_iter(text);
        if matches.is_empty() {
            return text.to_string();
        }
        let chars = &matches[0].text;
        let mut out = String::with_capacity(text.len());
        let mut last = 0;
        for found in &matches {
            out.extend(&chars[last..found.start()]);
            out.push_str(replacement);
            last = found.end();
        }
        out.extend(&chars[last..]);
        out
    }

    /// `re.split(pattern, text)` for a pattern without groups.
    pub fn split(&self, text: &str) -> Vec<String> {
        let chars: Vec<char> = text.chars().collect();
        let mut parts = Vec::new();
        let mut last = 0;
        for found in self.find_iter(text) {
            parts.push(chars[last..found.start()].iter().collect());
            last = found.end();
        }
        parts.push(chars[last..].iter().collect());
        parts
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
