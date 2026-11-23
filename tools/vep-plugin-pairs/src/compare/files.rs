// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The file pairs of a comparison over two directories.

use super::text::path_str;
use std::collections::BTreeMap;
use std::fs;

/// The endings of an output file's name, tried in this order for the pairing name.
const SUFFIXES: [&str; 6] = [".vep.vcf.gz", ".vep.vcf", ".vcf.gz", ".vcf", ".txt", ".gz"];
/// The endings a directory's entry is taken for an output by, in the order the entries of
/// each are listed.
const PATTERNS: [&str; 5] = [".txt", ".vep.vcf", ".vep.vcf.gz", ".vcf", ".vcf.gz"];

/// The name two outputs of one input are paired under: the file's name without the first
/// of `SUFFIXES` it ends in and without a leading `perl_` or `rust_`.
pub fn canonical_stem(name: &str) -> &str {
    let name = SUFFIXES
        .iter()
        .find_map(|suffix| name.strip_suffix(suffix))
        .unwrap_or(name);
    name.strip_prefix("perl_")
        .or_else(|| name.strip_prefix("rust_"))
        .unwrap_or(name)
}

/// The outputs of a directory by pairing name. Entries are listed per ending of `PATTERNS`
/// in the directory's own order, and of two entries with one pairing name the one listed
/// last stands. A path that is no readable directory has none.
fn outputs(root: &str) -> BTreeMap<String, String> {
    let root = path_str(root);
    let names: Vec<String> = fs::read_dir(&root)
        .map(|entries| {
            entries
                .filter_map(|entry| entry.ok()?.file_name().into_string().ok())
                .collect()
        })
        .unwrap_or_default();
    let mut found = BTreeMap::new();
    for pattern in PATTERNS {
        for name in names.iter().filter(|name| name.ends_with(pattern)) {
            let path = match root.as_str() {
                "." => name.clone(),
                "/" | "//" => format!("{root}{name}"),
                _ => format!("{root}/{name}"),
            };
            found.insert(canonical_stem(name).to_string(), path);
        }
    }
    found
}

/// The file pairs of two directories: each pairing name both hold, in the names' order,
/// with the Ensembl VEP side's path, the vep-rs side's path and the name.
pub fn matching_pairs(perl_dir: &str, rust_dir: &str) -> Vec<(String, String, String)> {
    let perl = outputs(perl_dir);
    let mut rust = outputs(rust_dir);
    perl.into_iter()
        .filter_map(|(name, perl_path)| {
            let rust_path = rust.remove(&name)?;
            Some((perl_path, rust_path, name))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_pairing_name_drops_one_ending_and_the_engine_prefix() {
        let cases = [
            ("perl_chr21.vep.vcf.gz", "chr21"),
            ("rust_chr21.vep.vcf", "chr21"),
            ("chr21.vcf.gz", "chr21"),
            ("chr21.vcf", "chr21"),
            ("perl_chr21.txt", "chr21"),
            ("chr21.txt.gz", "chr21.txt"),
            ("perl_rust_x.txt", "rust_x"),
            ("rust_perl_x.txt", "perl_x"),
            ("perl_.txt", ""),
            (".txt", ""),
            ("plain", "plain"),
            ("x.vep.txt", "x.vep"),
            ("Perl_x.TXT", "Perl_x.TXT"),
        ];
        for (name, expected) in cases {
            assert_eq!(canonical_stem(name), expected, "{name}");
        }
    }
}
