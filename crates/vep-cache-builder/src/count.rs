// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A count and its noun as a log line words them: a count of one takes the
//! noun in the singular (`1 transcript`), every other count, zero included,
//! the plural (`0 transcripts`, `2 transcripts`).

use std::fmt;

/// A count beside its noun; see [`counted`] and [`counted_as`].
pub struct Counted<N> {
    count: N,
    singular: &'static str,
    /// The plural where it is not the singular with an `s`.
    plural: Option<&'static str>,
}

/// `count` and `noun`, displayed as `1 <noun>` or `<count> <noun>s`.
pub fn counted<N>(count: N, noun: &'static str) -> Counted<N> {
    Counted {
        count,
        singular: noun,
        plural: None,
    }
}

/// `count` and a noun whose plural is not the singular with an `s`,
/// displayed as `1 <singular>` or `<count> <plural>`.
pub fn counted_as<N>(count: N, singular: &'static str, plural: &'static str) -> Counted<N> {
    Counted {
        count,
        singular,
        plural: Some(plural),
    }
}

/// `one` for a count of one, else `more`: the verb or pronoun that agrees
/// with a counted noun.
pub fn one_or<N: PartialEq + From<u8>>(
    count: N,
    one: &'static str,
    more: &'static str,
) -> &'static str {
    if count == N::from(1) {
        one
    } else {
        more
    }
}

impl<N: fmt::Display + PartialEq + From<u8>> fmt::Display for Counted<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let count = &self.count;
        match self.plural {
            _ if *count == N::from(1) => write!(f, "{count} {}", self.singular),
            Some(plural) => write!(f, "{count} {plural}"),
            None => write!(f, "{count} {}s", self.singular),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_count_of_one_takes_the_singular_and_every_other_count_the_plural() {
        let shown = |count: usize| counted(count, "region file").to_string();
        assert_eq!(shown(0), "0 region files");
        assert_eq!(shown(1), "1 region file");
        assert_eq!(shown(2), "2 region files");
        assert_eq!(shown(11), "11 region files");
        assert_eq!(counted(1u64, "byte").to_string(), "1 byte");
        assert_eq!(counted(21u64, "byte").to_string(), "21 bytes");
    }

    #[test]
    fn a_noun_with_a_plural_of_its_own_is_shown_in_it() {
        let shown = |count: usize| counted_as(count, "matrix", "matrices").to_string();
        assert_eq!(shown(0), "0 matrices");
        assert_eq!(shown(1), "1 matrix");
        assert_eq!(shown(2), "2 matrices");
    }

    #[test]
    fn the_word_beside_a_count_agrees_with_it() {
        assert_eq!(one_or(0usize, "names", "name"), "name");
        assert_eq!(one_or(1usize, "names", "name"), "names");
        assert_eq!(one_or(2usize, "names", "name"), "name");
        assert_eq!(one_or(1u64, "was", "were"), "was");
    }
}
