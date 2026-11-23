// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One match of a pattern, held apart from the text it was found in.

/// One group of a match: its text and where it stands in the searched text, in bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group {
    pub text: String,
    pub start: usize,
    pub end: usize,
}

/// The groups of one match, the whole match first; a group that took no part is `None`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Match {
    pub groups: Vec<Option<Group>>,
    /// The named groups with their indexes.
    pub names: Vec<(String, usize)>,
}

impl Match {
    /// The match `captures` of `pattern`.
    pub fn of(pattern: &regex::Regex, captures: &regex::Captures<'_>) -> Match {
        let groups = captures
            .iter()
            .map(|group| {
                group.map(|found| Group {
                    text: found.as_str().to_string(),
                    start: found.start(),
                    end: found.end(),
                })
            })
            .collect();
        Match {
            groups,
            names: names(pattern.capture_names()),
        }
    }

    /// The match `captures` of a `pattern` with a lookaround or a backreference.
    pub fn of_fancy(
        pattern: &fancy_regex::Regex,
        captures: &fancy_regex::Captures<'_, str>,
    ) -> Match {
        let groups = captures
            .iter()
            .map(|group| {
                group.map(|found| Group {
                    text: found.as_str().to_string(),
                    start: found.start(),
                    end: found.end(),
                })
            })
            .collect();
        Match {
            groups,
            names: names(pattern.capture_names()),
        }
    }

    /// The text of group `index`, 0 for the whole match; `None` for a group that took no
    /// part or that the pattern lacks.
    pub fn group(&self, index: usize) -> Option<&str> {
        self.groups
            .get(index)
            .and_then(|group| group.as_ref())
            .map(|group| group.text.as_str())
    }

    /// The text of the group called `name`.
    pub fn named(&self, name: &str) -> Option<&str> {
        self.names
            .iter()
            .find(|(held, _)| held == name)
            .and_then(|(_, index)| self.group(*index))
    }
}

fn names<'a>(names: impl Iterator<Item = Option<&'a str>>) -> Vec<(String, usize)> {
    names
        .enumerate()
        .filter_map(|(index, name)| name.map(|name| (name.to_string(), index)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_match_keeps_its_groups_by_index_and_by_name() {
        let pattern = regex::Regex::new(r"(?P<first>a+)(b)?(?P<last>c)").unwrap();
        let found = Match::of(&pattern, &pattern.captures("xaac").unwrap());
        assert_eq!(found.group(0), Some("aac"));
        assert_eq!(found.named("first"), Some("aa"));
        assert_eq!(found.group(2), None);
        assert_eq!(found.named("last"), Some("c"));
        assert_eq!(
            found.groups[0]
                .as_ref()
                .map(|group| (group.start, group.end)),
            Some((1, 4))
        );
        assert_eq!(found.named("absent"), None);
    }

    #[test]
    fn a_match_of_a_pattern_with_a_lookbehind_reads_the_same() {
        let pattern = fancy_regex::Regex::new(r"(?<=\d)([+-])(\d+)").unwrap();
        let found = Match::of_fancy(&pattern, &pattern.captures("12+34").unwrap().unwrap());
        assert_eq!(found.group(0), Some("+34"));
        assert_eq!(found.group(1), Some("+"));
        assert_eq!(found.group(2), Some("34"));
    }
}
