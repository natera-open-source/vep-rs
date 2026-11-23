// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `info.txt` in, `info.json` out: the lines as Ensembl VEP's own reader takes them,
//! the two values read off the cache directory's path, and the contig order.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::json::{PerlStr, Pretty};
use crate::number::numify;

/// A byte string of `info.txt` or of a path as Perl characters in UTF-8.
fn characters(bytes: &[u8]) -> Vec<u8> {
    PerlStr::bytes(bytes).characters().into_owned()
}

/// The `info.json` document while the run builds it.
#[derive(Debug, Default)]
pub struct Info {
    /// Every key but `source_versions`.
    pub top: BTreeMap<Vec<u8>, Pretty>,
    /// The `source_<key>` lines, by key without the prefix.
    pub source_versions: BTreeMap<Vec<u8>, Vec<u8>>,
}

impl Info {
    /// The lines of `info.txt` as Ensembl VEP's own reader takes them: a line that
    /// starts with `#` or has no tab is skipped; `source_<key>` is a source version;
    /// `variation_cols` is a list; `cell_types` and `regulatory` are dropped; any other
    /// key is kept unless its value is `-`. A later line replaces an earlier one.
    pub fn from_lines(text: &[u8]) -> Info {
        let mut info = Info::default();
        let mut lines: Vec<&[u8]> = text.split(|byte| *byte == b'\n').collect();
        // A final newline ends the last line; it does not start an empty one.
        if lines.last() == Some(&&b""[..]) {
            lines.pop();
        }
        for line in lines {
            if line.first() == Some(&b'#') {
                continue;
            }
            let Some(tab) = line.iter().position(|byte| *byte == b'\t') else {
                continue;
            };
            let (key, value) = (&line[..tab], &line[tab + 1..]);
            if key == b"cell_types" {
                continue;
            }
            if let Some(source) = key.strip_prefix(b"source_") {
                info.source_versions
                    .insert(characters(source), characters(value));
            } else if key == b"variation_cols" {
                let columns = value
                    .split(|byte| *byte == b',')
                    .filter(|column| !column.is_empty())
                    .map(|column| Pretty::String(characters(column)))
                    .collect();
                info.top.insert(characters(key), Pretty::Array(columns));
            } else if key == b"regulatory" {
                continue;
            } else if value != b"-" {
                info.top
                    .insert(characters(key), Pretty::String(characters(value)));
            }
        }
        info
    }

    /// The value of a top-level key that is a string.
    pub fn string(&self, key: &[u8]) -> Option<&[u8]> {
        match self.top.get(key) {
            Some(Pretty::String(value)) => Some(value),
            _ => None,
        }
    }

    /// Whether `$info->{source_versions}{$tool}` is true: present and neither empty
    /// nor the string `0`.
    pub fn has_source_version(&self, tool: &[u8]) -> bool {
        self.source_versions
            .get(tool)
            .is_some_and(|version| !(version.is_empty() || version == b"0"))
    }

    /// The whole document.
    pub fn document(&self) -> Pretty {
        let mut members = self.top.clone();
        members.insert(
            b"source_versions".to_vec(),
            Pretty::Object(
                self.source_versions
                    .iter()
                    .map(|(key, value)| (key.clone(), Pretty::String(value.clone())))
                    .collect(),
            ),
        );
        Pretty::Object(members)
    }
}

/// The last component of a cache directory path when it is `<digits>` or
/// `<digits>_<assembly>`, and the path before it: what the two Perl patterns
/// `/(\d+)(?:_[^/]+)?/?$` and `(?:^|/)([^/]+)/\d+(?:_[^/]+)?/?$` read. One trailing `/`
/// is allowed, and one newline at the very end as Perl's `$` allows it.
fn version_component(path: &[u8]) -> Option<(&[u8], &[u8])> {
    let path = path.strip_suffix(b"\n").unwrap_or(path);
    let path = path.strip_suffix(b"/").unwrap_or(path);
    let slash = path.iter().rposition(|byte| *byte == b'/')?;
    let (parent, component) = (&path[..slash], &path[slash + 1..]);
    let digits = component
        .iter()
        .position(|byte| !byte.is_ascii_digit())
        .unwrap_or(component.len());
    if digits == 0 {
        return None;
    }
    match &component[digits..] {
        [] => {}
        // `_` and at least one more character; the component holds no `/`, and a
        // newline inside it is a character like another.
        [b'_', rest @ ..] if !rest.is_empty() => {}
        _ => return None,
    }
    Some((parent, &component[..digits]))
}

/// The cache version: the digits of the path's last component, as a number.
pub fn cache_version(cache_dir: &[u8]) -> Option<Pretty> {
    let (_, digits) = version_component(cache_dir)?;
    Some(match numify(digits).plus(crate::number::Num::Int(0)) {
        crate::number::Num::Int(value) => Pretty::Integer(value),
        other => Pretty::Double(other.as_f64()),
    })
}

/// The transcript set: `refseq` or `merged` when the species directory (the component
/// before the version) ends `_refseq` or `_merged`, else `ensembl`; `None` for a path
/// without a species directory.
pub fn source_type(cache_dir: &[u8]) -> Option<&'static [u8]> {
    let (parent, _) = version_component(cache_dir)?;
    let species = match parent.iter().rposition(|byte| *byte == b'/') {
        Some(slash) => &parent[slash + 1..],
        None => parent,
    };
    if species.is_empty() {
        return None;
    }
    // Perl's `$` also matches before a newline that ends the name.
    let name = species.strip_suffix(b"\n").unwrap_or(species);
    Some(if name.ends_with(b"_refseq") {
        b"refseq"
    } else if name.ends_with(b"_merged") {
        b"merged"
    } else {
        b"ensembl"
    })
}

/// The contig order: names of digits alone first, by number; then `X`,
/// `Y` and `MT`; then every other name by its bytes. The sort is stable, as Perl's is.
pub fn contig_order(contigs: &[Vec<u8>]) -> Vec<Vec<u8>> {
    fn tier(contig: &[u8]) -> u8 {
        if !contig.is_empty() && contig.iter().all(u8::is_ascii_digit) {
            0
        } else if contig == b"X" {
            1
        } else if contig == b"Y" {
            2
        } else if contig == b"MT" {
            3
        } else {
            4
        }
    }
    let mut ordered = contigs.to_vec();
    ordered.sort_by(|left, right| {
        let (left_tier, right_tier) = (tier(left), tier(right));
        if left_tier != right_tier {
            return left_tier.cmp(&right_tier);
        }
        if left_tier == 0 {
            numify(left)
                .compare(numify(right))
                .unwrap_or(Ordering::Equal)
        } else {
            left.cmp(right)
        }
    });
    ordered
}

#[cfg(test)]
#[path = "info_tests.rs"]
mod tests;
