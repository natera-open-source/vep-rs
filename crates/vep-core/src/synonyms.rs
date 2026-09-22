// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Chromosome synonym resolver.
//!
//! Handles mapping between different chromosome naming conventions
//! (e.g., "chr21" <-> "21", "chrM" <-> "MT").

use std::collections::HashMap;

/// Chromosome synonym resolver.
#[derive(Debug, Clone, Default)]
pub struct ChromosomeSynonyms {
    /// Map from synonym to canonical name.
    synonyms: HashMap<String, String>,
}

impl ChromosomeSynonyms {
    pub fn new() -> Self {
        Self::default()
    }

    /// Load default synonyms (chr prefix mapping).
    pub fn with_defaults() -> Self {
        let mut s = Self::new();
        for chr in (1..=22)
            .map(|n| n.to_string())
            .chain(["X", "Y", "MT", "M"].iter().map(|c| c.to_string()))
        {
            s.add_synonym(&format!("chr{}", chr), &chr);
            s.add_synonym(&chr, &format!("chr{}", chr));
        }
        s.add_synonym("chrM", "MT");
        s.add_synonym("M", "MT");
        s
    }

    /// Load synonyms from a file (tab-delimited: name\tsynonym).
    pub fn load_from_file(path: &std::path::Path) -> std::io::Result<Self> {
        let mut s = Self::new();
        let content = std::fs::read_to_string(path)?;
        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                s.add_synonym(parts[0], parts[1]);
                s.add_synonym(parts[1], parts[0]);
            }
        }
        Ok(s)
    }

    /// Add a synonym mapping.
    pub fn add_synonym(&mut self, from: &str, to: &str) {
        self.synonyms.insert(from.to_string(), to.to_string());
    }

    /// Resolve a chromosome name to its canonical form.
    /// Returns the input unchanged if no synonym found.
    pub fn resolve<'a>(&'a self, chr: &'a str) -> &'a str {
        self.synonyms.get(chr).map(|s| s.as_str()).unwrap_or(chr)
    }

    /// Check if a chromosome matches another considering synonyms.
    pub fn matches(&self, a: &str, b: &str) -> bool {
        a == b || self.resolve(a) == b || a == self.resolve(b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_default_synonyms() {
        let s = ChromosomeSynonyms::with_defaults();

        assert_eq!(s.resolve("chr21"), "21");
        assert_eq!(s.resolve("21"), "chr21");
        assert_eq!(s.resolve("chrX"), "X");
        assert_eq!(s.resolve("X"), "chrX");
        assert_eq!(s.resolve("chrM"), "MT");
        assert_eq!(s.resolve("M"), "MT");
    }

    #[test]
    fn test_matches() {
        let s = ChromosomeSynonyms::with_defaults();

        assert!(s.matches("chr21", "21"));
        assert!(s.matches("21", "chr21"));
        assert!(s.matches("21", "21"));
        assert!(s.matches("chrX", "X"));
        assert!(!s.matches("chr21", "22"));
    }

    #[test]
    fn test_resolve_unknown() {
        let s = ChromosomeSynonyms::with_defaults();
        assert_eq!(s.resolve("UNKNOWN"), "UNKNOWN");
    }

    #[test]
    fn test_load_from_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("synonyms.txt");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "# comment line").unwrap();
            writeln!(f, "NC_000021.9\t21").unwrap();
            writeln!(f, "NC_000023.11\tX").unwrap();
        }

        let s = ChromosomeSynonyms::load_from_file(&path).unwrap();
        assert_eq!(s.resolve("NC_000021.9"), "21");
        assert_eq!(s.resolve("21"), "NC_000021.9");
        assert_eq!(s.resolve("NC_000023.11"), "X");
    }

    #[test]
    fn test_add_synonym() {
        let mut s = ChromosomeSynonyms::new();
        s.add_synonym("foo", "bar");
        assert_eq!(s.resolve("foo"), "bar");
        assert_eq!(s.resolve("baz"), "baz");
    }

    #[test]
    fn test_empty_synonyms() {
        let s = ChromosomeSynonyms::new();
        assert_eq!(s.resolve("21"), "21");
        assert!(s.matches("21", "21"));
        assert!(!s.matches("21", "chr21"));
    }
}
