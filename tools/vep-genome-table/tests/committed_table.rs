// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The committed table under `tests/genomes`: it holds one row for every ready-made cache
//! of the pinned releases, and every row is exact or names why it is not. The variable
//! `GENOME_TABLE_DIR` points the three tests at another table directory.

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use vep_genome_table::{
    check, outcome_kind, OutcomeKind, Row, BACTERIA_FILE, GENOMES_FILE, OUTCOMES_FILE, SETS_FILE,
};

/// The caches of `pub/release-116/variation/indexed_vep_cache/` on the Ensembl site.
const VERTEBRATES: u64 = 524;
/// The caches of `pub/metazoa/release-63/variation/indexed_vep_cache/` on the Ensembl
/// Genomes site.
const METAZOA: u64 = 385;
/// The caches of `pub/plants/release-63/variation/indexed_vep_cache/` on the Ensembl
/// Genomes site.
const PLANTS: u64 = 267;
/// The caches of `pub/fungi/release-63/variation/indexed_vep_cache/` on the Ensembl
/// Genomes site, its collection folders left out.
const FUNGI: u64 = 67;
/// The caches of `pub/protists/release-63/variation/indexed_vep_cache/` on the Ensembl
/// Genomes site, its collection folders left out.
const PROTISTS: u64 = 43;
/// The assemblies with a cache under
/// `pub/rapid-release/species/Homo_sapiens/<accession>/ensembl/variation/2022_10/indexed_vep_cache/`
/// on the Ensembl site.
const PANGENOME: u64 = 89;
/// The folders `bacteria_<N>_collection` of
/// `pub/bacteria/release-63/variation/indexed_vep_cache/` on the Ensembl Genomes site.
const BACTERIAL_COLLECTIONS: usize = 129;
/// The caches in those bacterial folders.
const BACTERIAL_CACHES: u64 = 31_289;
/// The folders `fungi_<name>_collection` of
/// `pub/fungi/release-63/variation/indexed_vep_cache/` on the Ensembl Genomes site.
const FUNGAL_COLLECTIONS: usize = 14;
/// The caches in those fungal folders.
const FUNGAL_COLLECTION_CACHES: u64 = 1_432;
/// The folders `protists_<name>_collection` of
/// `pub/protists/release-63/variation/indexed_vep_cache/` on the Ensembl Genomes site.
const PROTIST_COLLECTIONS: usize = 13;
/// The caches in those protist folders.
const PROTIST_COLLECTION_CACHES: u64 = 193;

/// The sets named whole: set, release, cache version, caches.
const NAMED: [(&str, &str, &str, u64); 6] = [
    ("vertebrates", "116", "116", VERTEBRATES),
    ("metazoa", "63", "63", METAZOA),
    ("plants", "63", "63", PLANTS),
    ("fungi", "63", "63", FUNGI),
    ("protists", "63", "63", PROTISTS),
    ("pangenome", "2022_10", "107", PANGENOME),
];

/// The collection families: the prefix of their sets, their release and cache version,
/// how many sets and how many caches.
const FAMILIES: [(&str, &str, usize, u64); 3] = [
    ("bacteria_", "63", BACTERIAL_COLLECTIONS, BACTERIAL_CACHES),
    ("fungi_", "63", FUNGAL_COLLECTIONS, FUNGAL_COLLECTION_CACHES),
    (
        "protists_",
        "63",
        PROTIST_COLLECTIONS,
        PROTIST_COLLECTION_CACHES,
    ),
];

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// How many lines of a failure a message prints before it counts the rest.
const SHOWN: usize = 20;

fn table_dir() -> PathBuf {
    match env::var_os("GENOME_TABLE_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => {
            let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/genomes");
            dir.canonicalize().unwrap_or(dir)
        }
    }
}

fn bytes(dir: &Path, name: &str) -> Vec<u8> {
    let path = dir.join(name);
    fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The lines of a table file after its header line, each with its line number. A file
/// that starts with the gzip magic is decompressed.
fn body(dir: &Path, name: &str) -> Vec<(usize, String)> {
    let raw = bytes(dir, name);
    let mut text = String::new();
    let read = if raw.starts_with(&GZIP_MAGIC) {
        MultiGzDecoder::new(raw.as_slice()).read_to_string(&mut text)
    } else {
        raw.as_slice().read_to_string(&mut text)
    };
    if let Err(error) = read {
        panic!("{}: {error}", dir.join(name).display());
    }
    text.lines()
        .enumerate()
        .skip(1)
        .map(|(index, line)| (index + 1, line.to_owned()))
        .collect()
}

/// Every line of the two row files as `<file>:<line>` and the row it parses to, or what
/// the line breaks.
fn rows(dir: &Path) -> Vec<(String, Result<Row, String>)> {
    [GENOMES_FILE, BACTERIA_FILE]
        .into_iter()
        .flat_map(|name| {
            body(dir, name)
                .into_iter()
                .map(move |(line, text)| (format!("{name}:{line}"), Row::parse(&text)))
        })
        .collect()
}

/// The sets file as set -> (release, cache version, caches).
fn sets(dir: &Path) -> BTreeMap<String, (String, String, u64)> {
    body(dir, SETS_FILE)
        .into_iter()
        .map(|(line, text)| {
            let fields: Vec<&str> = text.split('\t').collect();
            let [set, release, cache_version, caches] = fields.as_slice() else {
                panic!("{SETS_FILE}:{line}: {} fields, not 4", fields.len());
            };
            let caches = caches
                .parse()
                .unwrap_or_else(|_| panic!("{SETS_FILE}:{line}: caches \"{caches}\""));
            (
                (*set).to_owned(),
                ((*release).to_owned(), (*cache_version).to_owned(), caches),
            )
        })
        .collect()
}

/// The first [`SHOWN`] of `lines`, one a line, and how many more there are.
fn shown(lines: &[String]) -> String {
    let mut text = lines
        .iter()
        .take(SHOWN)
        .map(|line| format!("  {line}"))
        .collect::<Vec<_>>()
        .join("\n");
    if lines.len() > SHOWN {
        text.push_str(&format!("\n  and {} more", lines.len() - SHOWN));
    }
    text
}

#[test]
fn the_committed_table_is_complete() {
    let dir = table_dir();
    if let Err(failures) = check(&dir, true) {
        let lines: Vec<String> = failures.iter().map(ToString::to_string).collect();
        panic!(
            "the table under {} does not hold exactly the caches of sets.tsv, each key once: \
             {} failures\n{}",
            dir.display(),
            lines.len(),
            shown(&lines)
        );
    }
}

#[test]
fn every_row_is_exact_or_names_its_exclusion() {
    let dir = table_dir();
    let declared: BTreeSet<String> = body(&dir, OUTCOMES_FILE)
        .into_iter()
        .filter_map(|(_, line)| line.split('\t').next().map(str::to_owned))
        .collect();
    let mut by_outcome: BTreeMap<String, usize> = BTreeMap::new();
    let mut failures = Vec::new();
    for (place, row) in rows(&dir) {
        let row = match row {
            Ok(row) => row,
            Err(what) => {
                failures.push(format!("{place}: {what}"));
                continue;
            }
        };
        *by_outcome.entry(row.outcome.clone()).or_default() += 1;
        match outcome_kind(&row.outcome) {
            Some(OutcomeKind::Exact) => {
                let digest =
                    row.sha256.len() == 64 && row.sha256.bytes().all(|b| b.is_ascii_hexdigit());
                let records = row.records.parse::<u64>().is_ok_and(|count| count > 0);
                if !digest || !records {
                    failures.push(format!(
                        "{place}: {} is exact with records {} and sha256 {}",
                        row.key(),
                        row.records,
                        row.sha256
                    ));
                }
            }
            Some(OutcomeKind::Differs | OutcomeKind::Excluded)
                if declared.contains(&row.outcome) => {}
            _ => failures.push(format!(
                "{place}: {} has the outcome {}, which {OUTCOMES_FILE} does not declare",
                row.key(),
                row.outcome
            )),
        }
    }
    // The rows that are not exact, counted under the name each one gives.
    let named: Vec<String> = by_outcome
        .iter()
        .filter(|(outcome, _)| outcome.as_str() != "exact")
        .map(|(outcome, count)| format!("{outcome}: {count} rows"))
        .collect();
    println!(
        "{} rows are exact; rows that are not, by name:\n{}",
        by_outcome.get("exact").copied().unwrap_or(0),
        shown(&named)
    );
    assert!(
        failures.is_empty(),
        "{} rows under {} are neither exact with a digest and records nor a declared \
         outcome\n{}\nrows that are not exact, by name:\n{}",
        failures.len(),
        dir.display(),
        shown(&failures),
        shown(&named)
    );
}

#[test]
fn the_sets_are_the_pinned_releases() {
    let dir = table_dir();
    let sets = sets(&dir);
    let mut failures = Vec::new();

    for (set, release, cache_version, caches) in NAMED {
        let expected = (release.to_owned(), cache_version.to_owned(), caches);
        match sets.get(set) {
            Some(listed) if *listed == expected => {}
            Some((release, cache_version, caches)) => failures.push(format!(
                "{SETS_FILE}: {set} reads release {release}, cache version {cache_version}, \
                 {caches} caches, not {} {} {}",
                expected.0, expected.1, expected.2
            )),
            None => failures.push(format!("{SETS_FILE}: no line for {set}")),
        }
    }
    let mut expected_sets = NAMED.len();
    for (prefix, version, collections, caches) in FAMILIES {
        expected_sets += collections;
        let family: Vec<_> = sets
            .iter()
            .filter(|(set, _)| set.starts_with(prefix) && set.ends_with("_collection"))
            .collect();
        let held: u64 = family.iter().map(|(_, listed)| listed.2).sum();
        if family.len() != collections || held != caches {
            failures.push(format!(
                "{SETS_FILE}: {} sets {prefix}*_collection with {held} caches, not \
                 {collections} with {caches}",
                family.len()
            ));
        }
        for (set, (release, cache_version, _)) in family {
            if release != version || cache_version != version {
                failures.push(format!(
                    "{SETS_FILE}: {set} reads release {release} and cache version \
                     {cache_version}, not {version}"
                ));
            }
        }
    }
    if sets.len() != expected_sets {
        failures.push(format!(
            "{SETS_FILE}: {} sets, not {expected_sets}",
            sets.len()
        ));
    }

    // A line that is no row is the finding of the test of the rows.
    for (place, row) in rows(&dir) {
        let Ok(row) = row else { continue };
        match sets.get(&row.set) {
            Some((_, cache_version, _)) if *cache_version == row.cache_version => {}
            Some((_, cache_version, _)) => failures.push(format!(
                "{place}: {} has cache version {}, and its set reads {cache_version}",
                row.key(),
                row.cache_version
            )),
            None => failures.push(format!(
                "{place}: {} names the set {}, which {SETS_FILE} lacks",
                row.key(),
                row.set
            )),
        }
    }

    if bytes(&dir, GENOMES_FILE).starts_with(&GZIP_MAGIC) {
        failures.push(format!("{GENOMES_FILE} is gzip, not plain text"));
    }
    if !bytes(&dir, BACTERIA_FILE).starts_with(&GZIP_MAGIC) {
        failures.push(format!("{BACTERIA_FILE} is not gzip"));
    }

    assert!(
        failures.is_empty(),
        "the table under {} differs from the pinned releases: {} failures\n{}",
        dir.display(),
        failures.len(),
        shown(&failures)
    );
}
